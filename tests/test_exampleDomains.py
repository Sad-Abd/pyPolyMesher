import numpy as np
import pytest

from pyPolyMesher.exampleDomains import (
    CookDomain,
    HornDomain,
    MbbDomain,
    MichellDomain,
    SuspensionDomain,
    WrenchDomain,
)


ALL_DOMAINS = [
    CookDomain,
    SuspensionDomain,
    MichellDomain,
    WrenchDomain,
    HornDomain,
    MbbDomain,
]

DOMAINS_WITH_BC = [
    CookDomain,
    SuspensionDomain,
    MichellDomain,
    WrenchDomain,
    MbbDomain,
]


@pytest.mark.parametrize("domain", ALL_DOMAINS, ids=lambda d: d.name)
def test_example_domain_distfnc_on_small_grid(domain):
    xmin, xmax, ymin, ymax = domain.BdBox
    x = np.linspace(xmin, xmax, 4)
    y = np.linspace(ymin, ymax, 4)
    X, Y = np.meshgrid(x, y)
    points = np.column_stack((X.ravel(), Y.ravel()))

    result = domain.DistFnc(points)

    assert result.shape[0] == points.shape[0]
    assert result.shape[1] >= 2
    assert np.isfinite(result).all()


@pytest.mark.parametrize("domain", DOMAINS_WITH_BC, ids=lambda d: d.name)
def test_example_domain_boundary_conditions_return_valid_support_and_load(domain):
    xmin, xmax, ymin, ymax = domain.BdBox
    x = np.linspace(xmin, xmax, 9)
    y = np.linspace(ymin, ymax, 9)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack((X.ravel(), Y.ravel()))

    supp, load = domain.BndryCnds(nodes)

    assert supp is not None
    assert load is not None
    assert isinstance(supp, np.ndarray)
    assert isinstance(load, np.ndarray)
    assert supp.ndim == 2 and supp.shape[1] == 3
    assert load.ndim == 2 and load.shape[1] == 3
    assert supp.shape[0] > 0
    assert load.shape[0] > 0


@pytest.mark.parametrize("domain", ALL_DOMAINS, ids=lambda d: d.name)
def test_example_domain_plot_smoke(domain):
    domain.Plot(n=20)


def test_example_domain_smoke_inventory_is_complete():
    assert [domain.name for domain in ALL_DOMAINS] == [
        "Cook Domain",
        "Suspension Domain",
        "Michell Domain",
        "Wrench Domain",
        "Horn Domain",
        "Mbb Domain",
    ]
