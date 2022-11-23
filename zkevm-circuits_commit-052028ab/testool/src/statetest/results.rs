use anyhow::Result;
use handlebars::Handlebars;
use prettytable::Row;
use prettytable::Table;
use serde::Deserialize;
use serde::Serialize;
use serde_json::json;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::Read;
use std::io::Write;
use std::path::PathBuf;
use std::str::FromStr;
use strum::IntoEnumIterator;
use strum_macros::{EnumIter, EnumString}; // 0.17.1

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, EnumIter, EnumString, Serialize, Deserialize)]
pub enum ResultLevel {
    Success,
    Ignored,
    Fail,
    Panic,
}

#[derive(Eq, PartialEq, Clone, Debug, Serialize, Deserialize)]
pub struct ResultInfo {
    pub level: ResultLevel,
    pub details: String,
    pub path: String,
}

impl ResultLevel {
    pub fn display_string(&self) -> String {
        use ResultLevel::*;
        match self {
            Panic => "💀PANIC",
            Fail => "🔴FAILD",
            Ignored => "🟠IGNOR",
            Success => "🟢SUCCS",
        }
        .to_string()
    }
}
pub struct DiffEntry {
    id: String,
    prev: Option<ResultInfo>,
    curr: Option<ResultInfo>,
}

pub struct Diffs {
    info: String,
    tests: Vec<DiffEntry>,
}

impl Diffs {
    pub fn gen_info(&self) -> Vec<String> {
        let mut stat: HashMap<ResultLevel, isize> = HashMap::new();
        let mut stat_news = 0isize;

        for t in &self.tests {
            if let Some(prev) = &t.prev {
                *stat.entry(prev.level).or_default() -= 1;
                *stat.entry(t.curr.as_ref().unwrap().level).or_default() += 1;
            } else {
                stat_news += 1;
            }
        }

        let mut buff = String::default();
        buff.push_str(&format!("new: {:+} ", stat_news));
        for (lvl, n) in stat {
            buff.push_str(&format!("/ {:?}: {:+} ", lvl, n));
        }

        let mut out = Vec::new();
        out.push(format!("DIFFS {}\n", self.info));
        out.push(buff);
        for t in &self.tests {
            if let Some(prev) = &t.prev {
                let curr = t.curr.as_ref().unwrap();
                out.push(format!(
                    "{} : {:?}({}) => {:?}({})\n",
                    t.id, prev.level, prev.details, curr.level, curr.details
                ));
            }
        }

        out
    }
}

pub struct Report {
    tests: HashMap<String, ResultInfo>,
    diffs: Diffs,
    by_folder: Table,
    by_result: Table,
}

impl Report {
    pub fn print_tty(&self) -> Result<()> {
        self.by_folder.print_tty(false)?;
        self.by_result.print_tty(false)?;
        println!("{:#?}", self.diffs.gen_info());
        Ok(())
    }
    pub fn gen_html(&self) -> Result<String> {
        let template = include_str!("report.handlebars");
        let reg = Handlebars::new();
        let mut by_folder = Vec::new();
        let mut by_result = Vec::new();

        self.by_folder.print_html(&mut by_folder)?;
        self.by_result.print_html(&mut by_result)?;

        let data = &json!({
                "by_folder": String::from_utf8(by_folder)?,
                "by_result" : String::from_utf8(by_result)? ,
                "diffs" : self.diffs.gen_info(),
                "all_results" : self.tests
        });

        let html = reg.render_template(template, data)?;
        Ok(html)
    }
}

#[derive(Default)]
pub struct Results {
    tests: HashMap<String, ResultInfo>,
    cache: Option<PathBuf>,
}

impl Results {
    pub fn from_file(path: PathBuf) -> Result<Self> {
        let mut file = std::fs::File::open(&path)?;
        let mut buf = String::new();
        file.read_to_string(&mut buf)?;
        let mut tests = HashMap::new();
        for line in buf.lines().filter(|l| l.len() > 1) {
            let mut split = line.splitn(4, ';');
            let level = ResultLevel::from_str(split.next().unwrap()).unwrap();
            let id = split.next().unwrap().to_string();
            let path = split.next().unwrap().to_string();
            let details = split.next().unwrap().to_string();
            tests.insert(
                id,
                ResultInfo {
                    level,
                    path,
                    details,
                },
            );
        }
        Ok(Self { cache: None, tests })
    }

    pub fn with_cache(path: PathBuf) -> Result<Self> {
        let tests = if path.exists() {
            Self::from_file(path.clone())?.tests
        } else {
            HashMap::new()
        };
        Ok(Self {
            tests,
            cache: Some(path),
        })
    }

    pub fn report(self, previous: Option<(String, Results)>) -> Report {
        // collect data
        let mut folders = HashSet::new();
        let mut results = HashSet::new();
        let mut count_by_folder_level: HashMap<String, usize> = HashMap::new();
        let mut count_by_result: HashMap<String, usize> = HashMap::new();

        let mut diffs = Diffs {
            info: String::default(),
            tests: Vec::new(),
        };
        let mut prev_results = None;
        if let Some((prev_info, p_results)) = previous {
            diffs.info = prev_info;
            prev_results = Some(p_results);
        }

        for (id, info) in &self.tests {
            let name = &info.path.rsplit_terminator('/').next().unwrap();
            let folder = &info.path[..info.path.len() - name.len() - 1];
            let result = format!("{:?}_{}", info.level, info.details);

            folders.insert(folder);
            results.insert(result.to_string());

            let key = format!("{}_{:?}", folder, info.level);
            *count_by_folder_level.entry(key).or_default() += 1;
            *count_by_result.entry(result).or_default() += 1;

            if let Some(prev_results) = &prev_results {
                if let Some(prev_info) = prev_results.tests.get(id) {
                    if info != prev_info {
                        diffs.tests.push(DiffEntry {
                            id: id.to_string(),
                            prev: Some(prev_info.clone()),
                            curr: Some(info.clone()),
                        });
                    }
                } else {
                    diffs.tests.push(DiffEntry {
                        id: id.to_string(),
                        prev: None,
                        curr: Some(info.clone()),
                    });
                }
            }
        }

        let mut folders: Vec<_> = folders.iter().collect();
        folders.sort();
        let mut results: Vec<_> = results.iter().collect();
        results.sort();

        // generate tables

        let mut by_folder = Table::new();
        let mut header = vec![String::from("path")];

        let levels: Vec<_> = ResultLevel::iter().collect();

        header.append(&mut levels.iter().map(|v| format!("{:?}", v)).collect());
        by_folder.add_row(Row::from_iter(header));

        let mut totals = vec![0usize; levels.len()];

        for folder in folders {
            let mut row = Vec::new();
            for i in 0..levels.len() {
                let key = format!("{}_{:?}", folder, levels[i]);
                let value = *count_by_folder_level.get(&key).unwrap_or(&0usize);
                row.push(value);
                totals[i] += value;
            }
            let sum: usize = row.iter().sum();
            let mut cells = vec![folder.to_string()];
            cells.append(
                &mut row
                    .iter()
                    .map(|n| format!("{} ({}%)", n, (100 * n) / sum))
                    .collect(),
            );
            by_folder.add_row(Row::from_iter(cells));
        }
        let sum: usize = totals.iter().sum();
        let mut cells = vec!["TOTAL".to_string()];
        if sum != 0 {
            cells.append(
                &mut totals
                    .iter()
                    .map(|n| format!("{} ({}%)", n, (100 * n) / sum))
                    .collect(),
            );
        }
        by_folder.add_row(Row::from_iter(cells));

        let mut by_result = Table::new();
        let mut info = Vec::new();
        for (result, count) in count_by_result {
            info.push((count, result));
        }

        info.sort_by(|a, b| b.0.cmp(&a.0));
        for entry in info.iter().take(25) {
            by_result.add_row(row![format!("{}", entry.0), entry.1]);
        }

        Report {
            tests: self.tests,
            by_folder,
            by_result,
            diffs,
        }
    }

    pub fn contains(&self, test: &str) -> bool {
        self.tests.contains_key(test)
    }

    #[allow(clippy::map_entry)]
    pub fn insert(&mut self, test_id: String, result: ResultInfo) -> Result<()> {
        if !self.tests.contains_key(&test_id) {
            log::info!(
                "{} {}/{} {}",
                result.level.display_string(),
                result.path,
                test_id,
                result.details
            );

            let entry = format!(
                "{:?};{};{};{}\n",
                result.level, test_id, result.path, result.details
            );
            if let Some(path) = &self.cache {
                std::fs::OpenOptions::new()
                    .read(true)
                    .write(true)
                    .create(true)
                    .append(true)
                    .open(path)?
                    .write_all(entry.as_bytes())?;
            }
            self.tests.insert(test_id, result);
        }

        Ok(())
    }
}
