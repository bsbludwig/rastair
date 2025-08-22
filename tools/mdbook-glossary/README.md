# Mdbook Glossary

Transforms a YAML file with glossary entries to a chapter in [`mdbook`](https://github.com/rust-lang/mdBook)
and lets you easily link to entries.

## Configuration

In you `book.toml`, add:

```toml
[preprocessor.glossary]
command = "mdbook-glossary"
glossary = "src/glossary.yaml"
```

## Glossary file

A YAML file that is a list of entries with these fields:

| Field name    | Description                                         | Example                                     |
| ------------- | --------------------------------------------------- | ------------------------------------------- |
| `key`         | Identifier                                          | "BAM"                                       |
| `group`       | Group this entry belongs to, optional               | "File Formats"                              |
| `long`        | Long title of the entry, optional                   | "Binary Alignment Map"                      |
| `short`       | Short version of entry, optional, defaults to `key` | (empty, defaults to "BAM")                  |
| `plural`      | Plural version of `short`, optional                 | (empty, so it defaults to "BAMs")           |
| `description` | Description of the entry, in Markdown               | "A binary version of the @SAM file format." |

The glossary chapter will be appended to the book.
It will contain the entries group by the `group` key (if present)
and sorted alphabetically.

## Referencing entries

Entries can be referenced in the book using the `@key` syntax.
Multiple alternatives are supported to have a good text flow:

| Syntax    | Description                                    | Example     |
| --------- | ---------------------------------------------- | ----------- |
| `@key`    | Creates a short link to the entry              | `[key](#)`  |
| `@Key`    | Uppercase first letter of the key              | `[Key](#)`  |
| `@key:pl` | Plural form of the key, if available           | `[keys](#)` |
| `@Key:pl` | Uppercase first letter, plural form of the key | `[Keys](#)` |
