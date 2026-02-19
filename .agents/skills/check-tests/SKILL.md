---
name: check-tests
description: Run all Rastair tests and report which ones are failing. Use when you need a quick overview of test status before starting work.
context: fork
---

# Check Test Status

Run the full Rastair test suite and report only the failing tests.

## Task

Run tests efficiently and report only failures:

```bash
# Get test summary and failures in one go
cargo test 2>&1 | grep -E "(test result:|test.*\.\.\. FAILED)"
```

Then format the output as:

```
Test Status: X passed, Y failed

Failing tests:
- test_name_1
- test_name_2
```

If all tests pass, report: "✅ All tests passing"

## Notes

- The grep filters output to just summary + failures (much less noise)
- Extract test names from lines like `test path::to::test_name ... FAILED`
- Keep it concise - no need to show compilation warnings or test details
