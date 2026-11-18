"""End-to-end rain-to-loss pipeline.

Usage:  python -m rain2loss.run --config config.yaml
Every step writes its result to the output folder; finished steps are skipped on the next run.
"""
import argparse
import json
import logging
import time
from contextlib import contextmanager
from pathlib import Path

import geopandas as gpd
import pandas as pd
import yaml

from rain2loss.calibration import calibrate
from rain2loss.evt import (annual_maxima, block_bootstrap_years, bootstrap_return_levels,
                           empirical_return_periods, fit_gev_lmom, return_level,
                           scale_largest_storm)
from rain2loss.hbv import elevation_zones, run_hbv
from rain2loss.lamah import load_attributes, load_hbv_inputs
from rain2loss.loss import (MAX_DAMAGE_CZ_TOTAL, average_annual_loss, loss_table,
                            sample_depth)
from rain2loss.metrics import evaluate

log = logging.getLogger("rain2loss")


@contextmanager
def step(name):
    """Log the start, end and duration of a pipeline step."""
    log.info("START %s", name)
    t0 = time.perf_counter()
    yield
    log.info("END   %s (%.1f s)", name, time.perf_counter() - t0)


def setup_logging(out_dir):
    """Write log messages to the screen and to run.log (force=True: spotpy sets up its own logger)."""
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.FileHandler(out_dir / "run.log", mode="w"),
                                  logging.StreamHandler()],
                        force=True)


def main(config_path):
    cfg = yaml.safe_load(Path(config_path).read_text())
    out = Path(cfg["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    setup_logging(out)
    log.info("config: %s", cfg)

    gid = cfg["gauge_id"]
    with step("load data"):
        attrs = load_attributes()
        data = load_hbv_inputs(gid, attrs).loc[: cfg["periods"]["end"]]
        to_m3s = attrs.loc[gid, "area_gov"] / 86.4               # mm/day -> m3/s
        zones = elevation_zones(attrs.loc[gid, "elev_mean"], attrs.loc[gid, "elev_std"],
                                cfg["hbv"]["n_zones"])
        cal = (data.index >= cfg["periods"]["calibration"][0]) & (data.index <= cfg["periods"]["calibration"][1])
        val = (data.index >= cfg["periods"]["validation"][0]) & (data.index <= cfg["periods"]["validation"][1])

    params_file = out / "hbv_params.json"
    with step("calibrate HBV"):
        if params_file.exists() and not cfg["hbv"].get("recalibrate", False):
            log.info("using existing %s", params_file)
            params = json.loads(params_file.read_text())
        else:
            params, _ = calibrate(cfg["hbv"]["algorithm"], data, data["qobs"], cal, zones,
                                  repetitions=cfg["hbv"]["repetitions"], seed=cfg["seed"])
            params_file.write_text(json.dumps(params, indent=2))

    with step("evaluate HBV"):
        q_sim = run_hbv(data, params, zones)["q_sim"]
        scores = pd.DataFrame([{"period": name, **evaluate(q_sim[mask], data["qobs"][mask])}
                               for name, mask in [("calibration", cal), ("validation", val)]])
        scores.to_csv(out / "hbv_scores.csv", index=False)
        log.info("scores:\n%s", scores.round(3).to_string(index=False))

    with step("flood frequency"):
        ams = annual_maxima(data["qobs"] * to_m3s)
        T = cfg["return_periods"]
        gev = fit_gev_lmom(ams["peak"])
        boot = bootstrap_return_levels(ams["peak"], T, seed=cfg["seed"])
        freq = pd.DataFrame({"T": T, "Q_m3s": return_level(gev, T),
                             "Q_p05": boot.quantile(0.05).to_numpy(),
                             "Q_p95": boot.quantile(0.95).to_numpy()})
        freq.to_csv(out / "flood_frequency.csv", index=False)

    with step("event set"):
        forcing = data[["prec", "temp", "pet"]]
        synthetic = block_bootstrap_years(forcing, cfg["event_set"]["n_years"], seed=cfg["seed"])
        if cfg["event_set"].get("storm_sigma"):
            synthetic, _ = scale_largest_storm(synthetic, sigma=cfg["event_set"]["storm_sigma"],
                                               seed=cfg["seed"])
        synthetic["q_m3s"] = run_hbv(synthetic, params, zones)["q_sim"].to_numpy() * to_m3s
        events = (synthetic[synthetic["syn_year"] > 1]           # year 1 is warm-up
                  .groupby("syn_year").agg(peak_m3s=("q_m3s", "max"),
                                           source_year=("source_year", "first")))
        events["return_period"] = empirical_return_periods(events["peak_m3s"])
        events.sort_values("peak_m3s", ascending=False).to_csv(out / "event_set.csv")

    with step("losses"):
        buildings = gpd.read_file(cfg["loss"]["buildings"])
        buildings["value_eur"] = buildings.to_crs(3035).area * MAX_DAMAGE_CZ_TOTAL
        points = buildings.copy()
        points["geometry"] = buildings.to_crs(3035).representative_point().to_crs(buildings.crs)
        depth = {rp: sample_depth(cfg["loss"]["hazard_pattern"].format(rp=rp), points)
                 for rp in cfg["loss"]["return_periods"]}
        table = loss_table(points, depth)
        table.to_csv(out / "loss_by_return_period.csv", index=False)
        aal = average_annual_loss(table["rp"], table["loss_eur"],
                                  protection_rp=cfg["loss"].get("protection_rp"))
        (out / "aal.json").write_text(json.dumps({"aal_eur": aal}, indent=2))
        log.info("AAL = %.0f EUR", aal)

    log.info("done; results in %s", out.resolve())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rain-to-loss pipeline")
    parser.add_argument("--config", default="config.yaml", help="path to the YAML configuration")
    main(parser.parse_args().config)