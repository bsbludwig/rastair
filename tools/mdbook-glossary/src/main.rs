use anyhow::{Context, Result, ensure};
use clap::Parser as _;
use mdbook::{
    BookItem,
    book::{Book, Chapter},
    preprocess::{CmdPreprocessor, Preprocessor, PreprocessorContext},
};
use smol_str::{SmolStr, format_smolstr};
use std::{
    collections::{BTreeMap, HashMap},
    fmt::Write,
    io::{self, BufReader, BufWriter},
    path::Path,
};
use tracing::{debug, info, warn};

#[derive(Debug, clap::Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Option<Command>,
}

#[derive(Debug, clap::Subcommand)]
enum Command {
    Supports { renderer: String },
    Process,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt().with_writer(std::io::stderr).init();

    let pre = GlossaryProcessor;

    let cli = Cli::parse();
    match cli.command {
        Some(Command::Supports { .. }) => {
            return Ok(()); // support all renderers
        }
        _ => {
            let (ctx, book) = CmdPreprocessor::parse_input(BufReader::new(io::stdin()))?;
            let processed_book = pre.run(&ctx, book)?;
            serde_json::to_writer(BufWriter::new(io::stdout()), &processed_book)?;
        }
    }

    Ok(())
}

struct GlossaryProcessor;

impl Preprocessor for GlossaryProcessor {
    fn name(&self) -> &str {
        "GlossaryProcessor"
    }

    fn run(&self, ctx: &PreprocessorContext, mut book: Book) -> Result<Book> {
        let path = ctx
            .config
            .get("preprocessor.glossary.glossary")
            .and_then(|path| path.as_str())
            .unwrap_or("src/glossary.yaml");
        let glossary_file = ctx.root.join(path);
        ensure!(
            glossary_file.exists(),
            "Glossary file `{glossary_file:?}` not found in the book root."
        );
        debug!(?glossary_file, "Using glossary file");

        let glossary = std::fs::read_to_string(&glossary_file)
            .with_context(|| format!("Failed to read glossary file `{glossary_file:?}`"))?;
        let glossary = serde_yaml::from_str::<Glossary>(&glossary)
            .with_context(|| format!("Failed to parse glossary file `{glossary_file:?}`"))?;

        book.sections.push(BookItem::Chapter(Chapter {
            name: "Glossary".into(),
            content: render(glossary.by_group_sorted())
                .context("failed to write glossary items to markdown")?,
            number: None,
            sub_items: vec![],
            path: Some("glossary.md".into()),
            source_path: Some("glossary.yaml".into()),
            parent_names: vec![],
        }));
        info!(entries = glossary.0.len(), "added glossary as chapter");

        let replacements =
            Replacements(glossary.0.iter().fold(HashMap::new(), |mut map, entry| {
                let key = &entry.key;
                let short = entry.short.as_ref().unwrap_or(key).to_owned();
                let key_capitalized = capitalize(&entry.key);
                let plural = entry.plural.clone().unwrap_or(format_smolstr!("{short}s"));

                map.insert(
                    entry.key.clone(),
                    Replacement { short: short.clone(), entry: entry.clone() },
                );
                map.insert(
                    key_capitalized.clone(),
                    Replacement { short: capitalize(&short), entry: entry.clone() },
                );
                map.insert(
                    format_smolstr!("{key}:pl"),
                    Replacement { short: plural.clone(), entry: entry.clone() },
                );
                map.insert(
                    format_smolstr!("{key_capitalized}:pl"),
                    Replacement { short: capitalize(&plural), entry: entry.clone() },
                );
                map
            }));

        // Go thru each book chapter and replace `@ref` references with link to glossary entries.
        // This includes the glossary as well!
        book.for_each_mut(|item| {
            let BookItem::Chapter(chapter) = item else {
                return;
            };
            match replace_references(chapter, path, &replacements) {
                Ok(content) => chapter.content = content,
                Err(error) => {
                    let error = format!("{error:#}");
                    warn!(error, chapter.name, "Failed to replace references in chapter")
                }
            }
        });

        Ok(book)
    }

    fn supports_renderer(&self, _renderer: &str) -> bool {
        true
    }
}

fn capitalize(s: &str) -> SmolStr {
    s.chars().take(1).flat_map(|c| c.to_uppercase()).chain(s.chars().skip(1)).collect()
}

#[derive(Debug, serde::Deserialize)]
struct Glossary(Vec<GlossaryEntry>);

impl Glossary {
    fn by_group_sorted(&self) -> BTreeMap<SmolStr, Vec<GlossaryEntry>> {
        self.0.iter().fold(BTreeMap::new(), |mut acc, entry| {
            let group = entry.group.clone().unwrap_or_default();
            acc.entry(group).or_default().push(entry.clone());
            acc
        })
    }
}

#[derive(Debug, Clone, serde::Deserialize)]
struct GlossaryEntry {
    key: SmolStr,
    group: Option<SmolStr>,
    long: Option<SmolStr>,
    short: Option<SmolStr>,
    plural: Option<SmolStr>,
    description: SmolStr,
}

fn render(glossary_groups: BTreeMap<SmolStr, Vec<GlossaryEntry>>) -> Result<String> {
    let mut res = String::new();

    writeln!(res, "# Glossary")?;

    for (group, entries) in glossary_groups {
        writeln!(res, "## {group}")?;

        for entry in entries {
            entry.to_markdown(&mut res)?;
        }
    }

    Ok(res)
}

impl GlossaryEntry {
    fn to_markdown(&self, mut f: impl Write) -> Result<()> {
        let GlossaryEntry { key, group: _, long, short, plural: _, description } = self;
        let short = short.as_ref().unwrap_or(key);

        write!(f, "### **{short}**")?;
        if let Some(long) = long {
            write!(f, " – {long}")?;
        }
        write!(f, "  {{#entry--{key}}}")?;
        writeln!(f)?;

        writeln!(f, "{description}")?;

        Ok(())
    }
}

struct Replacements(HashMap<SmolStr, Replacement>);
struct Replacement {
    short: SmolStr,
    entry: GlossaryEntry,
}

fn replace_references(
    chapter: &Chapter,
    glossary_base_path: &str,
    replacements: &Replacements,
) -> Result<String> {
    use regex::Regex;

    let relative_glossary_path = {
        let base = glossary_base_path.replace(".yaml", ".md");
        if let Some(chapter_path) = chapter.path.as_ref() {
            // find relative path from chapter_path to glossary_base_path
            // e.g. given `src/foo/bar.md` and `src/glossary.md`, it is `../glossary.md`
            let chapter_dir = chapter_path.parent().unwrap_or_else(|| Path::new(""));
            let glossary_path = Path::new(&base);
            let glossary_path = glossary_path.strip_prefix("src").unwrap_or(glossary_path);

            if let Some(relative_path) = pathdiff::diff_paths(glossary_path, chapter_dir)
                && let Some(relative_path) = relative_path.to_str()
            {
                relative_path.into()
            } else {
                base
            }
        } else {
            base
        }
    };
    let entry_link = |key: &str| format!("{relative_glossary_path}#entry--{key}");

    let re = Regex::new(r"@([a-zA-Z0-9_-]+(?::pl)?)")?;
    let result = re.replace_all(&chapter.content, |caps: &regex::Captures| {
        let pattern = &caps[1];
        let key = pattern.strip_suffix(":pl").unwrap_or(pattern);

        if replacements.0.contains_key(pattern) {
            let x = &replacements.0[pattern];
            let display_text = &x.short;
            let hover_text = if let Some(long) = x.entry.long.as_ref() {
                format!(" \"{}\"", long.replace("\"", "\\\""))
            } else {
                "".into()
            };
            let link = entry_link(key);
            format!("[{display_text}]({link}{hover_text})")
        } else {
            caps[0].into() // Return original if not found in replacements
        }
    });

    Ok(result.into_owned())
}
