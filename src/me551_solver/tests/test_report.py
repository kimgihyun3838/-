"""Tests for ReportEngine – render_markdown and save_markdown.

These tests verify the report engine against SolverResult objects.
They do NOT depend on any solver implementation, so all tests should pass now.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from me551_solver.core.base import SolverResult
from me551_solver.core.report import ReportEngine


# ---------------------------------------------------------------------------
# render_markdown – structure tests
# ---------------------------------------------------------------------------

class TestRenderMarkdownStructure:
    """Verify section headers and ordering."""

    def test_title_is_problem_type(self, minimal_result):
        """Markdown title must be '# <problem_type>'."""
        md = ReportEngine.render_markdown(minimal_result)
        assert md.startswith("# Test Problem"), (
            "First line must be '# Test Problem'"
        )

    def test_given_section_present(self, minimal_result):
        """## Given section must appear in output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "## Given" in md

    def test_derivation_steps_section_present(self, minimal_result):
        """## Derivation Steps section must appear in output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "## Derivation Steps" in md

    def test_final_answer_section_present(self, minimal_result):
        """## Final Answer section must appear in output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "## Final Answer" in md

    def test_section_order(self, minimal_result):
        """Sections must appear in order: Given → Steps → Final Answer."""
        md = ReportEngine.render_markdown(minimal_result)
        idx_given = md.index("## Given")
        idx_steps = md.index("## Derivation Steps")
        idx_answer = md.index("## Final Answer")
        assert idx_given < idx_steps < idx_answer, (
            "Sections must be in order: Given, Steps, Final Answer"
        )

    def test_sanity_check_section_present_when_nonempty(self, minimal_result):
        """## Sanity Check section must appear when sanity_check is non-empty."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "## Sanity Check" in md

    def test_sanity_check_section_absent_when_empty(self):
        """## Sanity Check section must be absent when sanity_check is empty."""
        result = SolverResult(
            problem_type="No Check",
            given={},
            steps=[],
            final_answer={},
            sanity_check="",
        )
        md = ReportEngine.render_markdown(result)
        assert "## Sanity Check" not in md


# ---------------------------------------------------------------------------
# render_markdown – content tests
# ---------------------------------------------------------------------------

class TestRenderMarkdownContent:
    """Verify that given values, steps, and answers appear correctly."""

    def test_given_keys_and_values_rendered(self, minimal_result):
        """Each given key-value pair must appear as a list item."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "**m**" in md
        assert "1.0" in md
        assert "**k**" in md
        assert "100.0" in md

    def test_step_descriptions_rendered(self, minimal_result):
        """Step descriptions must appear as sub-headings."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "EOM" in md
        assert "omega_n" in md

    def test_step_expressions_in_math_blocks(self, minimal_result):
        """Step expressions must be wrapped in $$ ... $$ blocks."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "$$" in md
        assert "m x'' + k x = 0" in md

    def test_step_numbering(self, minimal_result):
        """Steps must be numbered starting from 1."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "Step 1:" in md
        assert "Step 2:" in md

    def test_final_answer_keys_rendered(self, minimal_result):
        """Final answer keys must appear bold in the output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "**omega_n**" in md

    def test_final_answer_values_rendered(self, minimal_result):
        """Final answer values must appear in the output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "10.0" in md

    def test_sanity_check_content_rendered(self, minimal_result):
        """Sanity check string must appear verbatim in output."""
        md = ReportEngine.render_markdown(minimal_result)
        assert "omega_n = sqrt(100/1) = 10 rad/s" in md

    def test_full_result_rendering(self, full_result):
        """Full FRF result must render all four steps and final answers."""
        md = ReportEngine.render_markdown(full_result)
        assert "Step 1:" in md
        assert "Step 4:" in md
        assert "G1_pole" in md
        assert "FRF Decomposition" in md


# ---------------------------------------------------------------------------
# render_markdown – edge cases
# ---------------------------------------------------------------------------

class TestRenderMarkdownEdgeCases:
    """Edge cases: empty dicts, no steps, special characters."""

    def test_empty_given(self):
        """Empty given dict must not crash."""
        result = SolverResult(
            problem_type="Empty Given",
            given={},
            steps=[("step", "x=0")],
            final_answer={"x": 0},
        )
        md = ReportEngine.render_markdown(result)
        assert "## Given" in md

    def test_empty_steps(self):
        """Empty steps list must not crash."""
        result = SolverResult(
            problem_type="No Steps",
            given={"a": 1},
            steps=[],
            final_answer={"a": 1},
        )
        md = ReportEngine.render_markdown(result)
        assert "## Derivation Steps" in md

    def test_empty_final_answer(self):
        """Empty final_answer dict must not crash."""
        result = SolverResult(
            problem_type="No Answer",
            given={"a": 1},
            steps=[("s", "e")],
            final_answer={},
        )
        md = ReportEngine.render_markdown(result)
        assert "## Final Answer" in md

    def test_returns_string(self, minimal_result):
        """render_markdown must always return a str."""
        result = ReportEngine.render_markdown(minimal_result)
        assert isinstance(result, str)

    def test_multiple_final_answer_entries(self):
        """All final answer entries must appear in output."""
        result = SolverResult(
            problem_type="Multi",
            given={},
            steps=[],
            final_answer={"alpha": 1.0, "beta": 2.0, "gamma": 3.0},
        )
        md = ReportEngine.render_markdown(result)
        for key in ("alpha", "beta", "gamma"):
            assert f"**{key}**" in md


# ---------------------------------------------------------------------------
# save_markdown
# ---------------------------------------------------------------------------

class TestSaveMarkdown:
    """Verify file I/O behaviour of save_markdown."""

    def test_creates_file(self, minimal_result):
        """save_markdown must create a .md file at the specified path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.md"
            ReportEngine.save_markdown(minimal_result, out)
            assert out.exists(), "Markdown file was not created"

    def test_file_content_matches_render(self, minimal_result):
        """File content must equal render_markdown output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.md"
            ReportEngine.save_markdown(minimal_result, out)
            expected = ReportEngine.render_markdown(minimal_result)
            actual = out.read_text(encoding="utf-8")
            assert actual == expected

    def test_creates_parent_directories(self, minimal_result):
        """save_markdown must create intermediate directories if they don't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "nested" / "deep" / "report.md"
            ReportEngine.save_markdown(minimal_result, out)
            assert out.exists()

    def test_accepts_string_path(self, minimal_result):
        """save_markdown must accept a plain str filepath."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = str(Path(tmpdir) / "report.md")
            ReportEngine.save_markdown(minimal_result, out)
            assert Path(out).exists()

    def test_file_encoding_utf8(self, minimal_result):
        """File must be written as UTF-8 (sanity check string may contain non-ASCII)."""
        result = SolverResult(
            problem_type="한국어 테스트",
            given={"m": 1.0},
            steps=[("EOM", "방정식")],
            final_answer={"ω_n": 10.0},
            sanity_check="확인 완료 ✓",
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report_utf8.md"
            ReportEngine.save_markdown(result, out)
            content = out.read_text(encoding="utf-8")
            assert "한국어 테스트" in content
            assert "확인 완료" in content

    def test_overwrites_existing_file(self, minimal_result, full_result):
        """Calling save_markdown twice must overwrite the first file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir) / "report.md"
            ReportEngine.save_markdown(minimal_result, out)
            ReportEngine.save_markdown(full_result, out)
            content = out.read_text(encoding="utf-8")
            assert "FRF Decomposition" in content
            assert "Test Problem" not in content
