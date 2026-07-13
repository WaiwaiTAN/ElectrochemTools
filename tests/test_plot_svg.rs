use electrochem_tools::plot::{PlotSeries, write_xy_svg};
use std::fs;

#[test]
fn svg_is_a_single_viewable_file_with_latex_labels_and_semantic_ids() {
    let root = std::env::temp_dir().join(format!(
        "electrochem_tools_latex_svg_{}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let path = root.join("plot.svg");

    write_xy_svg(
        &path,
        "A & <B> \"C\" 'D'",
        "Current density / mA cm^-2",
        "Frequency / Hz",
        &[PlotSeries {
            name: "A&B",
            points: vec![(1.0, 2.0), (2.0, 3.0)],
            color: "#123456",
        }],
    )
    .unwrap();

    let files = fs::read_dir(&root)
        .unwrap()
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert_eq!(files.len(), 1);
    assert_eq!(files[0].file_name(), "plot.svg");

    let svg = fs::read_to_string(path).unwrap();
    assert!(svg.contains(r#"<text id="title""#));
    assert!(svg.contains("A &amp; &lt;B&gt; &quot;C&quot; &apos;D&apos;"));
    assert!(svg.contains(r#"id="xlabel""#));
    assert!(svg.contains(r"$j$ / \unit{\milli\ampere\per\square\centi\metre}"));
    assert!(svg.contains(r#"id="ylabel""#));
    assert!(svg.contains(r"$f$ / \unit{\hertz}"));
    assert!(svg.contains(r#"id="legend-a-b""#));
    assert!(svg.contains("A&amp;B"));
    assert!(!svg.contains(">Current density / mA cm^-2<"));
    assert!(!svg.contains(">Frequency / Hz<"));
}
