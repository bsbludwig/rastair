# Quality control report

Rastair comes with a helper script, written in [R](https://cran.r-project.org), which generates a comprehense QC report for your library.

```admonish info
The QC tool requires a working installation of R with [RMarkdown](https://cran.r-project.org/web/packages/rmarkdown/index.html), [data.table](https://r-datatable.com) and [ggplot2](https://ggplot2.tidyverse.org) libraries.
```
## Generating QC reports
To generate the report, you need to first generate per-read bed output as described in the [examples section](src/examples.md#3-report-methylation-per-read).

```admonish tip
You can speed this up by e.g. restricting to a smaller chromosome with e.g. `-l chr17` as an additional argument to `rastair per-read`.
```

Once you have your per-read output, you generate the html report with

```bash
mkdir -p test_qc
mbias.R --output-prefix test_qc test_per-read.bed.gz
```

This will produce a file called `mbias.html` in the `test_qc` directory.

## Elements of the QC report

### M-bias

### V-bias