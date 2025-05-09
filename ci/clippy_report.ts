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

interface ClippyMessage {
  code: {
    code: string;
    explanation: string;
  };
  level: string;
  message: string;
  spans: Array<{
    file_name: string;
    line_start: number;
    line_end: number;
    column_start: number;
    column_end: number;
    is_primary: boolean;
    text: Array<{
      text: string;
      highlight_start: number;
      highlight_end: number;
    }>;
  }>;
  rendered: string;
}

interface ClippyOutput {
  reason: string;
  message: ClippyMessage;
}

async function runClippy(): Promise<ClippyOutput[]> {
  $.log("Running Clippy analysis...");

  const lines = await captureOutput([
    "cargo",
    "clippy",
    "--message-format=json",
    "--workspace",
    "--all-targets",
  ]);

  return lines
    .map((line) => {
      try {
        return JSON.parse(line) as ClippyOutput;
      } catch {
        return null;
      }
    })
    .filter(
      (item): item is ClippyOutput =>
        item !== null &&
        item.reason === "compiler-message" &&
        item.message?.level !== null &&
        ["warning", "error"].includes(item.message.level),
    );
}

async function main() {
  try {
    // Run Clippy and collect diagnostics
    const diagnostics = await runClippy();
    // Filter to separate errors and warnings
    const errors = diagnostics.filter((diag) => diag.message.level === "error");
    const warnings = diagnostics.filter((diag) =>
      diag.message.level === "warning"
    );
    $.log(
      `Processing ${diagnostics.length} Clippy diagnostics (${errors.length} errors, ${warnings.length} warnings)...`,
    );

    // Create report data
    const reportData: BitbucketReport = {
      title: "Rust Clippy Lint Results",
      details:
        `Clippy found ${diagnostics.length} warnings or errors in the codebase.`,
      report_type: "TEST",
      reporter: "clippy",
      result: diagnostics.length > 0 ? "FAILED" : "PASSED",
      data: [
        {
          title: "Total Issues",
          type: "NUMBER",
          value: diagnostics.length,
        },
      ],
    };

    // Submit main report
    await submitReport("clippy-report", reportData);

    // Process each diagnostic
    for (const output of diagnostics) {
      const message = output.message;
      const primarySpan = message.spans.find((span) => span.is_primary);

      if (primarySpan) {
        const fileName = primarySpan.file_name;
        const lineStart = primarySpan.line_start;
        const severity = message.level === "error" ? "HIGH" : "MEDIUM";
        const annotationId = createAnnotationId(
          "clippy",
          `${fileName}-${lineStart}`,
        );

        const annotationData: BitbucketAnnotation = {
          external_id: annotationId,
          title: `Clippy ${message.level}`,
          annotation_type: "CODE_SMELL",
          summary: message.rendered.trim(),
          severity,
          path: fileName,
          line: lineStart,
        };

        $.log(`Processing issue in ${fileName}:${lineStart}...`);
        await submitAnnotation("clippy-report", annotationId, annotationData);
      }
    }

    $.log("Clippy report processing complete!");

    // Exit with error only if compilation errors were found
    if (errors.length > 0) {
      $.logError(`${errors.length} errors found.`);
      Deno.exit(1);
    } else if (warnings.length > 0) {
      $.logWarn(`${warnings.length} warnings found, but proceeding.`);
    }
  } catch (error) {
    $.logError("Error processing Clippy report:", (error as Error).message);
    Deno.exit(1);
  }
}

if (import.meta.main) {
  await main();
}
