"""
Utility functions for single-cell analysis.

This module provides helper functions and utilities.
"""

import logging
import yaml
import os
from typing import Dict, Any

logger = logging.getLogger(__name__)


def setup_logging(log_file: str = "logs/analysis.log", level: str = "INFO") -> None:
    """
    Setup logging configuration.

    Parameters
    ----------
    log_file : str
        Path to log file
    level : str
        Logging level
    """
    # Create logs directory if it doesn't exist
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )


def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from YAML file.

    Parameters
    ----------
    config_path : str
        Path to configuration file

    Returns
    -------
    dict
        Configuration dictionary
    """
    logger.info(f"Loading configuration from {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    logger.info(f"Configuration loaded successfully")
    return config


def save_config(config: Dict[str, Any], config_path: str) -> None:
    """
    Save configuration to YAML file.

    Parameters
    ----------
    config : dict
        Configuration dictionary
    config_path : str
        Path to save configuration
    """
    os.makedirs(os.path.dirname(config_path), exist_ok=True)

    with open(config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    logger.info(f"Configuration saved to {config_path}")


def create_directories(output_dir: str) -> None:
    """
    Create necessary output directories.

    Parameters
    ----------
    output_dir : str
        Base output directory
    """
    subdirs = ["figures", "tables", "reports", "logs"]

    for subdir in subdirs:
        path = os.path.join(output_dir, subdir)
        os.makedirs(path, exist_ok=True)
        logger.info(f"Created directory: {path}")
