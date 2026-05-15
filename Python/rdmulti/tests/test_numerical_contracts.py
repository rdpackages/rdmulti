import contextlib
import io
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rdmulti-matplotlib"))

import numpy as np
import pandas as pd
import pytest

from rdmulti import rdmc, rdmcplot, rdms


DATA_DIR = Path(__file__).resolve().parents[2]


def quiet_call(fn, *args, **kwargs):
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        return fn(*args, **kwargs)


def unbox(value):
    return value[0] if isinstance(value, tuple) else value


def multicutoff_data():
    data_path = DATA_DIR / "simdata_multic.csv"
    if not data_path.exists():
        pytest.skip("repository illustration data not available")
    data = pd.read_csv(data_path)
    return data["y"], data["x"], data["c"]


def test_rdmc_matches_illustration_baseline():
    y, x, c = multicutoff_data()
    result = quiet_call(rdmc, y, x, c)

    np.testing.assert_allclose(
        unbox(result.Coefs).to_numpy(dtype=float),
        np.array([[484.83090323, 297.98107945, 398.91490949, 436.40049184]]),
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        unbox(result.Pv).to_numpy(dtype=float),
        np.array([[7.84407236e-48, 8.54527540e-16, 0.0, 7.41987105e-04]]),
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        unbox(result.H).to_numpy(dtype=float),
        np.array(
            [
                [14.6619828, 11.9522678, np.nan, 13.6840172],
                [14.6619828, 11.9522678, np.nan, 13.6840172],
            ]
        ),
        rtol=1e-8,
        atol=1e-8,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        unbox(result.Nh).to_numpy(dtype=float),
        np.array([[149.0, 120.0, 269.0, 273.0], [140.0, 126.0, 266.0, 277.0]]),
    )
    np.testing.assert_allclose(
        unbox(result.W),
        np.array([[0.54018692, 0.45981308]]),
        rtol=1e-8,
        atol=1e-8,
    )


def test_rdms_matches_cumulative_cutoff_baseline():
    data_path = DATA_DIR / "simdata_cumul.csv"
    if not data_path.exists():
        pytest.skip("repository illustration data not available")
    data = pd.read_csv(data_path)
    cvec = np.array([data["c"][0], data["c"][1]])
    result = quiet_call(rdms, data["y"], data["x"], cvec)

    np.testing.assert_allclose(
        result["Coefs"].to_numpy(dtype=float),
        np.array([[395.4918261, 342.87249624, np.nan]]),
        rtol=1e-8,
        atol=1e-8,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result["Pv"].to_numpy(dtype=float),
        np.array([[1.59995659e-145, 3.42385007e-120, np.nan]]),
        rtol=1e-8,
        atol=1e-8,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result["H"].to_numpy(dtype=float),
        np.array([[15.10927597, 12.22129226, np.nan], [15.10927597, 12.22129226, np.nan]]),
        rtol=1e-8,
        atol=1e-8,
        equal_nan=True,
    )
    np.testing.assert_allclose(
        result["Nh"].to_numpy(dtype=float),
        np.array([[140.0, 142.0, np.nan], [146.0, 123.0, np.nan]]),
        equal_nan=True,
    )


def test_rdmcplot_illustration_data_shape_is_stable():
    y, x, c = multicutoff_data()
    result = quiet_call(rdmcplot, y, x, c, nodraw=True)

    assert result["clist"].tolist() == [33.0, 66.0]
    assert result["Xmean"].shape == (2000, 2)
    assert result["Ymean"].shape == (2000, 2)
    assert result["X0"].shape == (2000, 2)
    assert result["X1"].shape == (2000, 2)
    np.testing.assert_allclose(
        result["Xmean"].head(5).to_numpy(dtype=float),
        np.array(
            [
                [0.53786151, 0.43492039],
                [1.76309871, 1.35887772],
                [2.78515738, 2.23753276],
                [3.85971366, 2.97432617],
                [4.78789688, 3.90511753],
            ]
        ),
        rtol=1e-6,
        atol=1e-6,
    )
