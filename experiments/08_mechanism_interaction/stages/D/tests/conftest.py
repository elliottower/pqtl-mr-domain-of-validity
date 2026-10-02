import pytest

from stage_d.join import derive, join_stages, load_stage_tables
from stage_d.sets import form_sets
from stage_d.synthetic import make_tables, write_stage_outputs


@pytest.fixture
def synthetic_tables():
    return make_tables(seed=3)


@pytest.fixture
def sealed_stages(tmp_path, synthetic_tables):
    root = tmp_path / "stages"
    expected = write_stage_outputs(root, *synthetic_tables)
    return root, expected


@pytest.fixture
def joined(sealed_stages):
    root, _ = sealed_stages
    return derive(join_stages(*load_stage_tables(root)))


@pytest.fixture
def analysis_sets(joined):
    df, _ = joined
    return form_sets(df)
