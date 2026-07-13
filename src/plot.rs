use anyhow::{Context, Result, bail};
use num_complex::Complex;
use std::fs;
use std::path::Path;

const WIDTH: f64 = 900.0;
const HEIGHT: f64 = 720.0;
const LEFT: f64 = 92.0;
const RIGHT: f64 = 40.0;
const TOP: f64 = 56.0;
const BOTTOM: f64 = 86.0;

pub struct PlotSeries<'a> {
    pub name: &'a str,
    pub points: Vec<(f64, f64)>,
    pub color: &'a str,
}

#[derive(Debug, Clone)]
struct Axis {
    min: f64,
    max: f64,
    ticks: Vec<Tick>,
}

#[derive(Debug, Clone)]
struct Tick {
    value: f64,
    label: String,
}

pub fn write_drt_gamma_svg(path: &Path, tau: &[f64], gamma: &[f64], title: &str) -> Result<()> {
    if tau.len() != gamma.len() || tau.is_empty() {
        bail!("cannot plot DRT gamma: tau/gamma length mismatch or empty data");
    }
    let points: Vec<(f64, f64)> = tau
        .iter()
        .zip(gamma)
        .filter_map(|(&tau, &gamma)| {
            (tau > 0.0 && tau.is_finite() && gamma.is_finite()).then_some((tau.log10(), gamma))
        })
        .collect();
    if points.is_empty() {
        bail!("cannot plot DRT gamma: no finite positive tau values");
    }
    let x_min = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
    let x_max = points
        .iter()
        .map(|(x, _)| *x)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_max_data = points
        .iter()
        .map(|(_, y)| y.max(0.0))
        .fold(0.0_f64, f64::max);
    let y_axis = nice_axis_from_zero(y_max_data, 5);
    let first_tick = x_min.ceil() as i32;
    let last_tick = x_max.floor() as i32;
    let x_axis = Axis {
        min: x_min,
        max: x_max.max(x_min + 1.0e-12),
        ticks: (first_tick..=last_tick)
            .map(|power| Tick {
                value: power as f64,
                label: format!("10^{power}"),
            })
            .collect(),
    };
    write_svg(
        path,
        title,
        "tau / s",
        "gamma",
        &[PlotSeries {
            name: "gamma",
            points,
            color: "#2563eb",
        }],
        x_axis,
        y_axis,
        false,
    )
}

pub fn write_nyquist_svg(
    path: &Path,
    title: &str,
    z_real: &[f64],
    z_imag: &[f64],
    z_fit: &[Complex<f64>],
) -> Result<()> {
    let experimental: Vec<(f64, f64)> = z_real
        .iter()
        .zip(z_imag)
        .filter_map(|(&re, &im)| (re.is_finite() && im.is_finite()).then_some((re, -im)))
        .collect();
    let fit: Vec<(f64, f64)> = z_fit
        .iter()
        .filter_map(|z| (z.re.is_finite() && z.im.is_finite()).then_some((z.re, -z.im)))
        .collect();
    if experimental.is_empty() && fit.is_empty() {
        bail!("cannot plot Nyquist data: no finite points");
    }
    let series = [
        PlotSeries {
            name: "experimental",
            points: experimental,
            color: "#111827",
        },
        PlotSeries {
            name: "fit",
            points: fit,
            color: "#dc2626",
        },
    ];
    let (x_axis, y_axis) = equal_ohm_axes(&series, 7);
    write_svg(
        path,
        title,
        "Z real / ohm",
        "-Z imag / ohm",
        &series,
        x_axis,
        y_axis,
        true,
    )
}

pub fn write_xy_svg(
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[PlotSeries<'_>],
) -> Result<()> {
    if series.is_empty() || series.iter().all(|s| s.points.is_empty()) {
        bail!("cannot plot empty data");
    }
    let (x_min, x_max, y_min, y_max) = bounds(series);
    let x_axis = nice_axis(x_min, x_max, 6, false);
    let y_axis = nice_axis(y_min, y_max, 6, false);
    write_svg(path, title, x_label, y_label, series, x_axis, y_axis, false)
}

#[allow(clippy::too_many_arguments)]
fn write_svg(
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[PlotSeries<'_>],
    x_axis: Axis,
    y_axis: Axis,
    square_plot: bool,
) -> Result<()> {
    let plot_size = if square_plot {
        (HEIGHT - TOP - BOTTOM).min(WIDTH - LEFT - RIGHT)
    } else {
        WIDTH - LEFT - RIGHT
    };
    let plot_w = plot_size;
    let plot_h = if square_plot {
        plot_size
    } else {
        HEIGHT - TOP - BOTTOM
    };
    let x_span = (x_axis.max - x_axis.min).abs().max(1.0e-12);
    let y_span = (y_axis.max - y_axis.min).abs().max(1.0e-12);
    let map = |x: f64, y: f64| {
        let px = LEFT + (x - x_axis.min) / x_span * plot_w;
        let py = TOP + (y_axis.max - y) / y_span * plot_h;
        (px, py)
    };

    let mut svg = String::new();
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">
<rect width="100%" height="100%" fill="white"/>
<text x="{cx}" y="30" text-anchor="middle" font-family="Arial" font-size="20">{title}</text>
<rect x="{LEFT}" y="{TOP}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#333" stroke-width="1"/>
"##,
        cx = WIDTH / 2.0,
        title = escape_xml(title)
    ));

    for tick in &x_axis.ticks {
        let (x, _) = map(tick.value, y_axis.min);
        svg.push_str(&format!(
            r##"<line x1="{x:.2}" y1="{TOP:.2}" x2="{x:.2}" y2="{y2:.2}" stroke="#e6e6e6"/>
<text x="{x:.2}" y="{ty:.2}" text-anchor="middle" font-family="Arial" font-size="11">{label}</text>
"##,
            y2 = TOP + plot_h,
            ty = TOP + plot_h + 20.0,
            label = escape_xml(&tick.label)
        ));
    }
    for tick in &y_axis.ticks {
        let (_, y) = map(x_axis.min, tick.value);
        svg.push_str(&format!(
            r##"<line x1="{LEFT:.2}" y1="{y:.2}" x2="{x2:.2}" y2="{y:.2}" stroke="#e6e6e6"/>
<text x="{tx:.2}" y="{yy:.2}" text-anchor="end" font-family="Arial" font-size="11">{label}</text>
"##,
            x2 = LEFT + plot_w,
            tx = LEFT - 8.0,
            yy = y + 4.0,
            label = escape_xml(&tick.label)
        ));
    }

    svg.push_str(&format!(
        r#"<text x="{cx}" y="{ly}" text-anchor="middle" font-family="Arial" font-size="14">{x_label}</text>
<text x="22" y="{cy}" text-anchor="middle" font-family="Arial" font-size="14" transform="rotate(-90 22 {cy})">{y_label}</text>
"#,
        cx = LEFT + plot_w / 2.0,
        ly = TOP + plot_h + 58.0,
        cy = TOP + plot_h / 2.0,
        x_label = escape_xml(x_label),
        y_label = escape_xml(y_label)
    ));

    for s in series {
        let points = s
            .points
            .iter()
            .filter(|(x, y)| x.is_finite() && y.is_finite())
            .map(|&(x, y)| {
                let (px, py) = map(x, y);
                format!("{px:.2},{py:.2}")
            })
            .collect::<Vec<_>>()
            .join(" ");
        svg.push_str(&format!(
            r#"<polyline points="{points}" fill="none" stroke="{color}" stroke-width="2"/>
"#,
            color = s.color
        ));
    }

    let mut legend_y = TOP + 20.0;
    for s in series {
        svg.push_str(&format!(
            r#"<line x1="{x1}" y1="{legend_y}" x2="{x2}" y2="{legend_y}" stroke="{color}" stroke-width="3"/>
<text x="{tx}" y="{ty}" font-family="Arial" font-size="12">{name}</text>
"#,
            x1 = LEFT + plot_w - 150.0,
            x2 = LEFT + plot_w - 120.0,
            tx = LEFT + plot_w - 112.0,
            ty = legend_y + 4.0,
            color = s.color,
            name = escape_xml(s.name)
        ));
        legend_y += 20.0;
    }

    svg.push_str("</svg>\n");
    fs::write(path, svg).with_context(|| format!("failed to write {}", path.display()))
}

fn equal_ohm_axes(series: &[PlotSeries<'_>], target_ticks: usize) -> (Axis, Axis) {
    let (mut x_min, mut x_max, mut y_min, mut y_max) = bounds(series);
    x_min = x_min.min(0.0);
    y_min = y_min.min(0.0);
    x_max = x_max.max(0.0);
    y_max = y_max.max(0.0);
    let span = (x_max - x_min).max(y_max - y_min).max(1.0e-12);
    let step = nice_step(span / target_ticks.max(1) as f64);
    let x_min_n = (x_min / step).floor() * step;
    let y_min_n = (y_min / step).floor() * step;
    let span_n = (span / step).ceil() * step;
    let x_axis = Axis {
        min: x_min_n,
        max: x_min_n + span_n,
        ticks: ticks_for_range(x_min_n, x_min_n + span_n, step),
    };
    let y_axis = Axis {
        min: y_min_n,
        max: y_min_n + span_n,
        ticks: ticks_for_range(y_min_n, y_min_n + span_n, step),
    };
    (x_axis, y_axis)
}

fn nice_axis(min: f64, max: f64, target_ticks: usize, force_zero_min: bool) -> Axis {
    let min = if force_zero_min { 0.0 } else { min };
    let span = (max - min).abs().max(1.0e-12);
    let step = nice_step(span / target_ticks.max(1) as f64);
    let axis_min = if force_zero_min {
        0.0
    } else {
        (min / step).floor() * step
    };
    let axis_max = (max / step).ceil() * step;
    Axis {
        min: axis_min,
        max: axis_max.max(axis_min + step),
        ticks: ticks_for_range(axis_min, axis_max.max(axis_min + step), step),
    }
}

fn nice_axis_from_zero(max: f64, target_ticks: usize) -> Axis {
    nice_axis(0.0, max.max(1.0e-12), target_ticks, true)
}

fn nice_step(raw_step: f64) -> f64 {
    if raw_step <= 0.0 || !raw_step.is_finite() {
        return 1.0;
    }
    let exponent = raw_step.log10().floor();
    let base = 10_f64.powf(exponent);
    let fraction = raw_step / base;
    let nice_fraction = if fraction <= 1.0 {
        1.0
    } else if fraction <= 2.0 {
        2.0
    } else if fraction <= 5.0 {
        5.0
    } else {
        10.0
    };
    nice_fraction * base
}

fn ticks_for_range(min: f64, max: f64, step: f64) -> Vec<Tick> {
    let mut ticks = Vec::new();
    let mut value = (min / step).ceil() * step;
    let limit = max + step * 1.0e-9;
    while value <= limit {
        ticks.push(Tick {
            value,
            label: format_tick(value),
        });
        value += step;
    }
    ticks
}

fn bounds(series: &[PlotSeries<'_>]) -> (f64, f64, f64, f64) {
    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    let mut y_min = f64::INFINITY;
    let mut y_max = f64::NEG_INFINITY;
    for s in series {
        for &(x, y) in &s.points {
            if x.is_finite() && y.is_finite() {
                x_min = x_min.min(x);
                x_max = x_max.max(x);
                y_min = y_min.min(y);
                y_max = y_max.max(y);
            }
        }
    }
    if x_min == x_max {
        x_min -= 0.5;
        x_max += 0.5;
    }
    if y_min == y_max {
        y_min -= 0.5;
        y_max += 0.5;
    }
    (x_min, x_max, y_min, y_max)
}

fn escape_xml(input: &str) -> String {
    input
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn format_tick(value: f64) -> String {
    let rounded = if value.abs() < 1.0e-12 { 0.0 } else { value };
    if (rounded - rounded.round()).abs() < 1.0e-9 {
        format!("{:.0}", rounded)
    } else if (rounded * 2.0 - (rounded * 2.0).round()).abs() < 1.0e-9 {
        format!("{rounded:.1}")
    } else if rounded.abs() >= 1.0e4 || (rounded != 0.0 && rounded.abs() < 1.0e-3) {
        format!("{rounded:.2e}")
    } else {
        format!("{rounded:.3}")
    }
}
