"""Pending review results must remain distinct from unsuccessful plans and errors."""

import json

from aquascope.gym import plans as gp


def _row(case_id, *, score, declined=False, correct=None, error=None, repeat=0):
    return gp.PlanResult(case_id=case_id, playbook="groundwater_decline", decline_expected=False,
                         agent="methodologist", model="example", provider="example", score=score,
                         declined=declined, decline_correct=correct, error=error, repeat=repeat)


def test_pending_rows_are_excluded_from_means_and_counted():
    rows = [_row("accepted", score=1.0, declined=True, correct=True),
            _row("rejected", score=0.0, declined=True, correct=False),
            _row("pending", score=None, declined=True),
            _row("zero-plan", score=0.0, correct=True),
            _row("partial-plan", score=0.5, correct=True),
            _row("error", score=0.0, error="offline")]
    summary = gp.summarize_plans(rows)[0]
    assert summary["score"] == 0.3  # includes the error, excludes pending
    assert summary["review_pending"] == 1
    assert summary["false_decline_rate"] == 0.25  # one of four judged, non-error cases
    assert summary["zero_plan_rate"] == 0.5  # one of two judged non-refusals
    assert summary["errors"] == 1
    report = gp.plan_leaderboard(rows)
    assert "review pending" in report and "zero plans" in report
    assert "excluded from score means" in report


def test_all_pending_has_no_score_and_old_json_still_loads(tmp_path):
    pending = _row("pending", score=None, declined=True)
    summary = gp.summarize_plans([pending])[0]
    assert summary["score"] is None and summary["by_playbook"] == {}
    assert summary["false_decline_rate"] is None and summary["review_pending"] == 1
    assert gp.PlanResult.from_dict(pending.to_dict()).score is None
    legacy = {"case_id": "legacy", "playbook": "x", "decline_expected": False,
              "agent": "tree", "model": None, "provider": None, "score": 0.9}
    path = tmp_path / "legacy.jsonl"
    path.write_text(json.dumps(legacy) + "\n", encoding="utf-8")
    assert gp.load_plan_results([path])[0].score == 0.9


def test_result_and_rescore_preserve_pending(monkeypatch):
    ref = gp.load_reference("gw_regional_potomac")
    candidate = gp.Candidate(declined=True, reason="ambiguous")
    result = {"score": None, "coverage_tools": None, "coverage_methods": None, "coverage_gates": None,
              "extraneous": None, "forbidden_used": 0, "decline_correct": None, "explain": ["review"]}
    monkeypatch.setattr(gp, "score_plan", lambda *_args, **_kwargs: result)
    row = gp._result(ref, "tree", gp._Config(), 0, candidate, {}, 0)
    assert row.score is None and row.decline_correct is None
    row.score = 0.0
    assert gp.rescore_plans([row], [ref])[0].score is None
