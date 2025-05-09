#!/usr/bin/env -S deno run --allow-all

import $ from "jsr:@david/dax";
import {
  BitbucketAnnotation,
  BitbucketReport,
  captureOutput,
  createAnnotationId,
  submitAnnotation,
  submitReport,
} from "./report_utils.ts";

interface TestResult {
  name: string;
  status: "passed" | "failed";
  output?: string;
  duration?: number;
}

async function runTests(): Promise<TestResult[]> {
  $.log("Running cargo tests...");
  let lines: string[] = [];
  try {
    lines = await captureOutput([
      "cargo",
      "test",
      "--all",
      "--no-fail-fast",
    ]);
  } catch (error) {
    // Attempt to extract output from the error if available
    if (error instanceof Error && "stdout" in (error as any)) {
      const output = new TextDecoder().decode((error as any).stdout);
      lines = output.split("\n").filter((line) => line.trim());
    } else {
      throw error;
    }
  }

  const results: TestResult[] = [];
  let currentTest: Partial<TestResult> & { output?: string[] } = {};

  for (const line of lines) {
    // Capture test results: "test test_name ... ok" or "test test_name ... FAILED"
    const resultLine = line.match(/^test\s+(\S+)\s+\.\.\.\s+(ok|FAILED)\s*$/);
    if (resultLine) {
      const [, testName, result] = resultLine;
      const status = result === "ok" ? "passed" : "failed";

      // If we already have captured output for this test, use it
      const output = currentTest.name === testName && currentTest.output?.length
        ? currentTest.output.join("\n")
        : undefined;

      results.push({ name: testName, status, ...(output && { output }) });
      currentTest = {};
      continue;
    }

    // Capture start of test with output: "test test_name ... "
    const testStart = line.match(/^test\s+(\S+)\s+\.\.\.\s*$/);
    if (testStart) {
      currentTest = { name: testStart[1], output: [] };
      continue;
    }

    // Capture test output
    if (currentTest.name && currentTest.output) {
      // Capture test failures with stack traces
      if (
        line.includes("thread '") &&
          (line.includes("' panicked at") || line.includes("' failed at")) ||
        line.startsWith("note: ") ||
        line.match(
          /^\s+left:|^\s+right:|^\s+at src\/.*:\d+$|^stack backtrace:/,
        ) ||
        // Also capture indented lines after a panic/failure for complete context
        (currentTest.output.length > 0 && line.startsWith("    "))
      ) {
        currentTest.output.push(line);
      }
    }
  }

  return results;
}

function parseTestOutput(output: string): { file?: string; line?: number } {
  // Match direct panic locations
  const panicMatch = output.match(
    /(?:panicked |failed )at .*?([^:\s]+\.rs):(\d+)(?::\d+)?/,
  );
  if (panicMatch) {
    return {
      file: panicMatch[1],
      line: parseInt(panicMatch[2], 10),
    };
  }

  // Match stack trace locations
  const stackMatch = output.match(/\bat (?:.*?[(`])?(src\/[^:]+):(\d+)/);
  if (stackMatch) {
    return {
      file: stackMatch[1],
      line: parseInt(stackMatch[2], 10),
    };
  }

  return {};
}

async function main() {
  try {
    // Run tests and parse results
    const results = await runTests();
    $.log(`Processing ${results.length} test results...`);

    // Count test results
    const passedTests = results.filter((r) => r.status === "passed").length;
    const failedTests = results.filter((r) => r.status === "failed").length;

    // Create report data
    const reportData: BitbucketReport = {
      title: "Rust Test Results",
      details: `${passedTests} tests passed, ${failedTests} tests failed.`,
      report_type: "TEST",
      reporter: "cargo-test",
      result: failedTests > 0 ? "FAILED" : "PASSED",
      data: [
        {
          title: "Total Tests",
          type: "NUMBER",
          value: results.length,
        },
        {
          title: "Passed Tests",
          type: "NUMBER",
          value: passedTests,
        },
        {
          title: "Failed Tests",
          type: "NUMBER",
          value: failedTests,
        },
      ],
    };

    // Submit main report
    await submitReport("test-report", reportData);

    // Submit annotations for failed tests
    for (const test of results) {
      if (test.status === "failed" && test.output) {
        const annotationId = createAnnotationId("test", test.name);
        const { file, line } = parseTestOutput(test.output);

        const annotationData: BitbucketAnnotation = {
          external_id: annotationId,
          title: `Test Failure: ${test.name}`,
          annotation_type: "BUG",
          summary: test.output.trim(),
          severity: "HIGH",
          ...(file && { path: file }),
          ...(line && { line }),
        };

        $.log(`Processing failure in ${test.name}...`);
        await submitAnnotation("test-report", annotationId, annotationData);
      }
    }

    $.log("Test report processing complete!");

    // Exit with error if any tests failed
    if (failedTests > 0) {
      $.logError(`${failedTests} tests failed.`);
      Deno.exit(1);
    }
  } catch (error) {
    $.logError("Error processing test report:", (error as Error).message);
    Deno.exit(1);
  }
}

if (import.meta.main) {
  await main();
}
