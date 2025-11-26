"""
Prepare scFEA inputs and invoke the scFEA CLI.

Usage:
    python src/scfea/run_scfea.py --config config/scfea_run_at2.yaml
    python src/scfea/run_scfea.py --config config/scfea_run_ipf.yaml

This assumes scFEA is installed in a conda env (default: py311) and that the
input counts CSV has already been exported from a Seurat object.
"""

import argparse
import datetime as dt
import os
import shutil
import subprocess
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scFEA with configurable inputs.")
    parser.add_argument("--config", type=str, help="Path to YAML config with parameters.")
    parser.add_argument("--counts_file", type=str, help="Path to the scRNA-seq counts matrix.")
    parser.add_argument("--input_dir", type=str, help="Directory passed to scFEA --input_dir.")
    parser.add_argument("--test_file_name", type=str, help="File name passed to scFEA --test_file.")
    parser.add_argument("--results_dir", type=str, help="Output directory for scFEA results.")
    parser.add_argument("--flux_output", type=str, help="Output CSV for predicted flux.")
    parser.add_argument("--balance_output", type=str, help="Output CSV for predicted balance.")
    parser.add_argument("--conda_env", type=str, help="Conda environment containing scFEA.")
    parser.add_argument("--scfea_root", type=str, help="Root directory of scFEA installation.")
    parser.add_argument("--scfea_data_dir", type=str, help="Directory passed to --data_dir.")
    parser.add_argument("--moduleGene_file", type=str, help="Module gene file name.")
    parser.add_argument("--stoichiometry_matrix", type=str, help="Stoichiometry matrix file name.")
    parser.add_argument("--cName_file", type=str, help="Compound name file name.")
    parser.add_argument("--sc_imputation", type=str, help="Whether to enable MAGIC imputation (True/False).")
    parser.add_argument("--train_epoch", type=int, help="Number of training epochs.")
    parser.add_argument("--log_file", type=str, help="Path to capture stdout/stderr from scFEA.")
    parser.add_argument("--log_dir", type=str, help="Directory for log files if --log_file not provided.")
    parser.add_argument("--use_gpu", type=str, help="Use GPU (True/False).")
    return parser.parse_args()


def load_config(path: str | None) -> dict:
    if not path:
        return {}
    cfg_path = Path(path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config file not found: {cfg_path}")
    with cfg_path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def coalesce(config: dict, args: argparse.Namespace, key: str, default=None):
    value = getattr(args, key, None)
    if value is not None:
        return value
    if key in config and config[key] is not None:
        return config[key]
    return default


def ensure_input_file(source: Path, destination_dir: Path, dest_name: str) -> Path:
    destination_dir.mkdir(parents=True, exist_ok=True)
    destination = destination_dir / dest_name
    if source.resolve() != destination.resolve():
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    return destination


def prepare_log_path(config: dict, args: argparse.Namespace) -> Path:
    log_file = coalesce(config, args, "log_file")
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        return log_path
    log_dir = Path(coalesce(config, args, "log_dir", "logs"))
    timestamp = dt.datetime.now().strftime("%Y%m%d%H%M%S")
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir / f"run-scfea-{timestamp}.log"


def ensure_patched_scfea_script(scfea_script: Path) -> Path:
    original_text = scfea_script.read_text(encoding="utf-8")
    patched_text = original_text.replace("detach().numpy()", "detach().cpu().numpy()")
    if patched_text == original_text:
        return scfea_script
    patched_path = scfea_script.with_name("scFEA_patched.py")
    patched_path.write_text(patched_text, encoding="utf-8")
    return patched_path


def normalize_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "1", "yes", "y"}:
            return True
        if lowered in {"false", "0", "no", "n"}:
            return False
    return bool(value)


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    counts_file_str = coalesce(config, args, "counts_file")
    if not counts_file_str:
        raise ValueError("counts_file must be provided via --counts_file or config.")
    counts_file = Path(counts_file_str).expanduser().resolve()
    input_dir = Path(coalesce(config, args, "input_dir", counts_file.parent)).expanduser().resolve()
    test_file_name = coalesce(config, args, "test_file_name", counts_file.name)
    ensured_counts = ensure_input_file(counts_file, input_dir, test_file_name)

    scfea_root = Path(coalesce(config, args, "scfea_root")).expanduser().resolve()
    scfea_script = scfea_root / "src" / "scFEA.py"
    if not scfea_script.exists():
        raise FileNotFoundError(f"scFEA script not found at {scfea_script}")
    scfea_exec_script = ensure_patched_scfea_script(scfea_script)
    scfea_data_dir = Path(coalesce(config, args, "scfea_data_dir", scfea_root / "data")).expanduser().resolve()

    results_dir = Path(coalesce(config, args, "results_dir", "results/scfea")).expanduser().resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    # Use absolute path to avoid issues with deep relative traversal
    results_dir_arg = str(results_dir)

    flux_output = Path(coalesce(config, args, "flux_output", results_dir / "adj_flux.csv")).expanduser().resolve()
    flux_output.parent.mkdir(parents=True, exist_ok=True)

    balance_output = Path(
        coalesce(config, args, "balance_output", results_dir / "adj_balance.csv")
    ).expanduser().resolve()
    balance_output.parent.mkdir(parents=True, exist_ok=True)
    module_gene_file = coalesce(config, args, "moduleGene_file", "module_gene_m168.csv")
    stoich_file = coalesce(config, args, "stoichiometry_matrix", "cmMat_c70_m168.csv")
    cname_file = coalesce(config, args, "cName_file", "cName_c70_m168.csv")
    batch_size = str(coalesce(config, args, "batch_size", 128))
    lr = str(coalesce(config, args, "lr", 0.01))
    sc_imputation_flag = normalize_bool(coalesce(config, args, "sc_imputation", True))
    train_epoch = str(coalesce(config, args, "train_epoch", 100))
    conda_env = coalesce(config, args, "conda_env", "py311")
    use_gpu = normalize_bool(coalesce(config, args, "use_gpu", False))

    log_path = prepare_log_path(config, args)

    cmd = [
        "conda",
        "run",
        "-n",
        conda_env,
        "python",
        str(scfea_exec_script),
        "--data_dir",
        str(scfea_data_dir),
        "--input_dir",
        str(input_dir),
        "--res_dir",
        results_dir_arg,
        "--test_file",
        test_file_name,
        "--moduleGene_file",
        module_gene_file,
        "--stoichiometry_matrix",
        stoich_file,
        "--cName_file",
        cname_file,
        "--sc_imputation",
        "True" if sc_imputation_flag else "False",
        "--output_flux_file",
        str(flux_output),
        "--output_balance_file",
        str(balance_output),
        "--train_epoch", train_epoch,
        "--batch_size", batch_size,
        "--lr", lr,
    ]

    log_path.write_text(f"Running command:\n{' '.join(cmd)}\n\n", encoding="utf-8")
    env = os.environ.copy()
    if not use_gpu:
        env["CUDA_VISIBLE_DEVICES"] = ""

    with log_path.open("a", encoding="utf-8") as log_handle:
        process = subprocess.run(
            cmd,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            cwd=str(scfea_root),
            env=env,
        )
    if process.returncode != 0:
        raise RuntimeError(f"scFEA execution failed (see log: {log_path})")

    if not flux_output.exists() or not balance_output.exists():
        raise FileNotFoundError(
          f"Expected scFEA outputs not found. Flux: {flux_output.exists()}, Balance: {balance_output.exists()}"
        )

    print(f"scFEA run completed. Flux: {flux_output}, Balance: {balance_output}, Log: {log_path}")


if __name__ == "__main__":
    main()
