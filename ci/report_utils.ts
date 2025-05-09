import $ from "jsr:@david/dax";

// Report types
export interface BitbucketReport {
  title: string;
  details: string;
  report_type: "SECURITY" | "COVERAGE" | "TEST" | "BUG";
  reporter: string;
  result: "PASSED" | "FAILED" | "PENDING";
  link?: string;
  logo_url?: string;
  data?: Array<{
    title: string;
    type?: "BOOLEAN" | "DATE" | "DURATION" | "LINK" | "NUMBER" | "PERCENTAGE" | "TEXT";
    value: number | boolean | string | { text: string; href: string };
  }>;
}

export interface BitbucketAnnotation {
  external_id: string;
  title?: string;
  annotation_type: "VULNERABILITY" | "CODE_SMELL" | "BUG";
  summary: string;
  details?: string;
  result?: "PASSED" | "FAILED" | "IGNORED" | "SKIPPED";
  severity: "HIGH" | "MEDIUM" | "LOW" | "CRITICAL";
  path?: string;
  line?: number;
  link?: string;
}

interface ReportContext {
  isLocal: boolean;
  repoOwner?: string;
  repoSlug?: string;
  commit?: string;
  outputDir: string;
}

// Get report context for either local or CI environment
function getReportContext(): ReportContext {
  const repoOwner = Deno.env.get("BITBUCKET_REPO_OWNER");
  const repoSlug = Deno.env.get("BITBUCKET_REPO_SLUG");
  const commit = Deno.env.get("BITBUCKET_COMMIT");
  const outputDir = Deno.env.get("LOCAL_REPORT_DIR") || "target/reports";

  const isLocal = !repoOwner || !repoSlug || !commit;

  if (isLocal) {
    $.log("Running in local development mode");
  } else {
    $.log("Running in Bitbucket Pipelines CI environment");
  }

  return { isLocal, repoOwner, repoSlug, commit, outputDir };
}

// Generate consistent IDs for annotations
export function createAnnotationId(prefix: string, identifier: string): string {
  return `${prefix}-${identifier}`.replace(/[^a-z0-9-]/gi, "-");
}

// Submit a report to either local file system or Bitbucket
export async function submitReport(
  reportKey: string,
  reportData: BitbucketReport,
): Promise<void> {
  const context = getReportContext();

  if (context.isLocal) {
    await writeLocalReport(reportKey, reportData, context);
  } else {
    await writeBitbucketReport(reportKey, reportData, context);
  }
}

// Submit an annotation to either local file system or Bitbucket
export async function submitAnnotation(
  reportKey: string,
  annotationId: string,
  annotationData: BitbucketAnnotation,
): Promise<void> {
  const context = getReportContext();

  if (context.isLocal) {
    await writeLocalAnnotation(
      reportKey,
      annotationId,
      annotationData,
      context,
    );
  } else {
    await writeBitbucketAnnotation(
      reportKey,
      annotationId,
      annotationData,
      context,
    );
  }
}

// Write local report file
async function writeLocalReport(
  reportKey: string,
  reportData: BitbucketReport,
  context: ReportContext,
): Promise<void> {
  await $.path(context.outputDir).mkdir({ recursive: true });
  const reportPath = $.path(context.outputDir).join(`${reportKey}.json`);

  $.log(`Writing report to ${reportPath}`);
  await Deno.writeTextFile(
    reportPath.toString(),
    JSON.stringify(reportData, null, 2),
  );
}

// Write local annotation file
async function writeLocalAnnotation(
  reportKey: string,
  annotationId: string,
  annotationData: BitbucketAnnotation,
  context: ReportContext,
): Promise<void> {
  const annotationsDir = $.path(context.outputDir).join(
    reportKey,
    "annotations",
  );
  await annotationsDir.mkdir({ recursive: true });

  const annotationPath = annotationsDir.join(`${annotationId}.json`);
  await Deno.writeTextFile(
    annotationPath.toString(),
    JSON.stringify(annotationData, null, 2),
  );
}

// Write report to Bitbucket using fetch with proxy
async function writeBitbucketReport(
  reportKey: string,
  reportData: BitbucketReport,
  context: ReportContext,
): Promise<void> {
  // When running in Pipelines, we need to use the proxy at localhost:29418
  // But use http instead of https as the proxy will handle the secure connection
  const url =
    `http://api.bitbucket.org/2.0/repositories/${context.repoOwner}/${context.repoSlug}/commit/${context.commit}/reports/${reportKey}`;

  try {
    // Set up the HTTP_PROXY environment variable for this request
    const originalHttpProxy = Deno.env.get("HTTP_PROXY");
    Deno.env.set("HTTP_PROXY", "http://localhost:29418");

    $.log(`Submitting report to ${url} through proxy http://localhost:29418`);

    const response = await fetch(url, {
      method: "PUT",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(reportData),
    });

    // Restore the original proxy setting if it existed
    if (originalHttpProxy) {
      Deno.env.set("HTTP_PROXY", originalHttpProxy);
    } else {
      Deno.env.delete("HTTP_PROXY");
    }

    if (!response.ok) {
      throw new Error(
        `HTTP error ${response.status}: ${await response.text()}`,
      );
    }

    $.log(`Report submitted successfully: ${response.status}`);
  } catch (error) {
    const errorMessage = error instanceof Error ? error.message : String(error);
    $.logError("Failed to submit report:", errorMessage);
  }
}

// Write annotation to Bitbucket using fetch with proxy
async function writeBitbucketAnnotation(
  reportKey: string,
  annotationId: string,
  annotationData: BitbucketAnnotation,
  context: ReportContext,
): Promise<void> {
  // When running in Pipelines, we need to use the proxy at localhost:29418
  const url =
    `http://api.bitbucket.org/2.0/repositories/${context.repoOwner}/${context.repoSlug}/commit/${context.commit}/reports/${reportKey}/annotations/${annotationId}`;

  try {
    // Set up the HTTP_PROXY environment variable for this request
    const originalHttpProxy = Deno.env.get("HTTP_PROXY");
    Deno.env.set("HTTP_PROXY", "http://localhost:29418");

    const location = annotationData.path && annotationData.line
      ? `${annotationData.path}:${annotationData.line}`
      : annotationId;

    $.log(
      `Submitting annotation for ${location} to ${url} through proxy http://localhost:29418`,
    );

    const response = await fetch(url, {
      method: "PUT",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(annotationData),
    });

    // Restore the original proxy setting if it existed
    if (originalHttpProxy) {
      Deno.env.set("HTTP_PROXY", originalHttpProxy);
    } else {
      Deno.env.delete("HTTP_PROXY");
    }

    if (!response.ok) {
      throw new Error(
        `HTTP error ${response.status}: ${await response.text()}`,
      );
    }

    $.log(
      `Annotation for ${location} submitted successfully: ${response.status}`,
    );
  } catch (error) {
    const location = annotationData.path && annotationData.line
      ? `${annotationData.path}:${annotationData.line}`
      : annotationId;
    const errorMessage = error instanceof Error ? error.message : String(error);
    $.logError(
      `Failed to submit annotation for ${location}:`,
      errorMessage,
    );
  }
}

// Run a process and capture its output lines while printing it in real-time
export async function captureOutput(
  command: string[],
  options: { cwd?: string } = {},
): Promise<string[]> {
  const cmdString = command.join(" ");
  $.log(`Running command: ${cmdString}`);

  const process = new Deno.Command(command[0], {
    args: command.slice(1),
    stdout: "piped",
    stderr: "piped",
    ...options,
  });

  const capturedLines: string[] = [];
  const child = process.spawn();

  // Handle output using events for simpler code
  const processStream = async (
    stream: ReadableStream<Uint8Array>,
    isError: boolean,
  ) => {
    const reader = stream.getReader();
    const decoder = new TextDecoder();

    try {
      while (true) {
        const { done, value } = await reader.read();
        if (done) break;

        const text = decoder.decode(value, { stream: true });
        const lines = text.split("\n");

        for (const line of lines) {
          if (line.trim()) {
            isError ? $.logError(line) : $.log(line);
            capturedLines.push(line);
          }
        }
      }
    } finally {
      reader.releaseLock();
    }
  };

  // Process both streams concurrently
  await Promise.all([
    processStream(child.stdout, false),
    processStream(child.stderr, true),
  ]);

  // Wait for process to complete
  await child.status;

  return capturedLines.filter((line) => line.trim());
}
