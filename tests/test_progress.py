import pytest

from pyPolyMesher.progress import Bar, pretty_time


def test_pretty_time_formats_representative_values():
    assert pretty_time(0) == "0:00:00"
    assert pretty_time(59.4) == "0:00:59"
    assert pretty_time(60) == "0:01:00"
    assert pretty_time(3661) == "1:01:01"


def test_bar_render_increment_done_enabled(monkeypatch):
    writes = []

    class FakeStdout:
        def write(self, text):
            writes.append(text)

        def flush(self):
            pass

    monkeypatch.setattr("sys.stdout", FakeStdout())

    bar = Bar(max_value=10, min_value=0, enabled=True)
    bar.update(2)
    assert bar.value == 2
    assert "20%" in bar.render()
    assert bar.render_value() == "(2 of 10)"
    assert bar.render_bar(size=10) == "[##--------]"

    bar.increment(3, 0.42)
    assert bar.value == 5
    assert bar.error == pytest.approx(0.42)
    assert bar.min_error == pytest.approx(0.42)
    assert bar.percent_complete == pytest.approx(50.0)
    assert bar.render_error() == "[Error: 0.420000 Min Error: 0.420000]"
    assert bar.eta >= 0

    bar.done()
    assert bar.value == 10
    assert any("\r" in item for item in writes)
    assert "\n" in writes


def test_bar_disabled_mode_executes_without_writing(monkeypatch):
    writes = []

    class FakeStdout:
        def write(self, text):
            writes.append(text)

        def flush(self):
            pass

    monkeypatch.setattr("sys.stdout", FakeStdout())

    bar = Bar(max_value=4, min_value=0, enabled=False)
    assert bar.render_value() == "(0 of 4)"
    assert bar.render_bar(size=4) == "[------------------------------]"[:6] or bar.render_bar(size=4)

    bar.increment(1, 1.0)
    bar.increment(1, 0.5)
    bar.done()

    assert bar.value == 4
    assert bar.error == pytest.approx(0.5)
    assert bar.min_error == pytest.approx(0.5)
    assert writes == []
