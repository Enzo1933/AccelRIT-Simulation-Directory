use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::{
    beam_and_tracker::{Beam, Tracker, beam_rigidity},
    magnet::MagnetGeometry,
};

// ============================================================
// App State
// ============================================================

pub struct QuadApp {
    // ── Beam ─────────────────────────────────────────────────
    drift_in: f64,
    energy_mev: f64,
    x0_mm: f64,
    xp0_mrad: f64,

    // ── Magnet Geometry ───────────────────────────────────────
    r_gap_mm: f64,   // bore radius
    l_mag_in: f64,   // magnet length
    w_pole_mm: f64,  // pole tip width
    l_iron_mm: f64,  // iron path length
    a_iron_mm2: f64, // iron cross section mm²
    mu_i: f64,       // initial permeability
    b_sat: f64,      // saturation field [T]
    gap: f64,        // Gap between magnets

    // ── Results ───────────────────────────────────────────────
    tracker: Option<Tracker>,
    mmf1: Option<f64>,
    mmf2: Option<f64>,
    g1: Option<f64>,
    g2: Option<f64>,

    // Flux results
    phi_total1: Option<f64>,
    phi_gap1: Option<f64>,
    phi_leak1: Option<f64>,
    phi_total2: Option<f64>,
    phi_gap2: Option<f64>,
    phi_leak2: Option<f64>,

    // Energy results
    e_total1: Option<f64>,
    e_gap1: Option<f64>,
    e_leak1: Option<f64>,
    e_total2: Option<f64>,
    e_gap2: Option<f64>,
    e_leak2: Option<f64>,

    status: String,

    // ── UI state ──────────────────────────────────────────────
    active_tab: Tab,
}

#[derive(PartialEq)]
enum Tab {
    Results,
    Fluxes,
    Energies,
}

impl Default for QuadApp {
    fn default() -> Self {
        Self {
            drift_in: 18.0,
            energy_mev: 1.0,
            x0_mm: 5.0,
            xp0_mrad: 30.0,

            r_gap_mm: 76.2,
            l_mag_in: 4.0,
            w_pole_mm: 60.0,
            l_iron_mm: 300.0,
            a_iron_mm2: 4000.0,
            mu_i: 2000.0,
            b_sat: 1.5,
            gap: 4.0,

            tracker: None,
            mmf1: None,
            mmf2: None,
            g1: None,
            g2: None,

            phi_total1: None,
            phi_gap1: None,
            phi_leak1: None,
            phi_total2: None,
            phi_gap2: None,
            phi_leak2: None,

            e_total1: None,
            e_gap1: None,
            e_leak1: None,
            e_total2: None,
            e_gap2: None,
            e_leak2: None,

            status: "Set parameters and press Run.".into(),
            active_tab: Tab::Results,
        }
    }
}

// ============================================================
// eframe App
// ============================================================

impl eframe::App for QuadApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // ── Top bar ──────────────────────────────────────────
        egui::TopBottomPanel::top("title").show(ctx, |ui| {
            ui.add_space(4.0);
            ui.heading("Quadrupole Triplet Design Tool — 1 MeV Proton Beam");
            ui.add_space(4.0);
        });

        // ── Left panel: inputs ────────────────────────────────
        egui::SidePanel::left("inputs")
            .min_width(300.0)
            .max_width(340.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    self.draw_inputs(ui);
                });
            });

        // ── Right panel: results tabs ─────────────────────────
        egui::SidePanel::right("results")
            .min_width(280.0)
            .max_width(320.0)
            .show(ctx, |ui| {
                self.draw_results(ui);
            });

        // ── Central panel: plot ───────────────────────────────
        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_plot(ui);
        });
    }
}

// ============================================================
// Panel drawing
// ============================================================

impl QuadApp {
    fn draw_inputs(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        // ── Beam section ──────────────────────────────────────
        ui.label(egui::RichText::new("Beam").strong().size(13.0));
        ui.separator();

        egui::Grid::new("beam_grid")
            .num_columns(2)
            .spacing([8.0, 6.0])
            .show(ui, |ui| {
                ui.label("Drift to target (in)");
                ui.add(egui::Slider::new(&mut self.drift_in, 6.0..=48.0).step_by(0.5));
                ui.end_row();

                ui.label("Energy (MeV)");
                ui.add(egui::Slider::new(&mut self.energy_mev, 0.1..=5.0).step_by(0.1));
                ui.end_row();

                ui.label("x₀ (mm)");
                ui.add(egui::Slider::new(&mut self.x0_mm, 0.5..=20.0).step_by(0.1));
                ui.end_row();

                ui.label("x′₀ (mrad)");
                ui.add(egui::Slider::new(&mut self.xp0_mrad, 1.0..=100.0).step_by(1.0));
                ui.end_row();
            });

        ui.add_space(12.0);

        // ── Magnet geometry section ───────────────────────────
        ui.label(egui::RichText::new("Magnet Geometry").strong().size(13.0));
        ui.separator();

        egui::Grid::new("magnet_grid")
            .num_columns(2)
            .spacing([8.0, 6.0])
            .show(ui, |ui| {
                ui.label("Bore radius (mm)");
                ui.add(egui::Slider::new(&mut self.r_gap_mm, 10.0..=150.0).step_by(1.0));
                ui.end_row();

                ui.label("Magnet length (in)");
                ui.add(egui::Slider::new(&mut self.l_mag_in, 1.0..=12.0).step_by(0.1));
                ui.end_row();

                ui.label("Pole width (mm)");
                ui.add(egui::Slider::new(&mut self.w_pole_mm, 10.0..=200.0).step_by(1.0));
                ui.end_row();

                ui.label("Iron path length (mm)");
                ui.add(egui::Slider::new(&mut self.l_iron_mm, 50.0..=1000.0).step_by(5.0));
                ui.end_row();

                ui.label("Iron area (mm²)");
                ui.add(egui::Slider::new(&mut self.a_iron_mm2, 100.0..=20000.0).step_by(100.0));
                ui.end_row();

                ui.label("μᵢ (initial permeability)");
                ui.add(egui::Slider::new(&mut self.mu_i, 100.0..=10000.0).step_by(100.0));
                ui.end_row();

                ui.label("B_sat (T)");
                ui.add(egui::Slider::new(&mut self.b_sat, 0.5..=2.5).step_by(0.05));
                ui.end_row();

                ui.label("Inter-magnet Gap(in)");
                ui.add(egui::Slider::new(&mut self.gap, 1.0..=100.0).step_by(0.1));
                ui.end_row();
            });

        ui.add_space(16.0);

        // ── Run button ────────────────────────────────────────
        if ui
            .add_sized(
                [280.0, 40.0],
                egui::Button::new(egui::RichText::new("Run optimizer").size(15.0)),
            )
            .clicked()
        {
            self.run();
        }

        ui.add_space(6.0);

        // Status line
        let status_color = if self.status.starts_with("Error") || self.status.starts_with("No") {
            egui::Color32::RED
        } else if self.tracker.is_some() {
            egui::Color32::GREEN
        } else {
            egui::Color32::GRAY
        };
        ui.label(egui::RichText::new(&self.status).color(status_color));

        // ── Export buttons ────────────────────────────────────
        if self.tracker.is_some() {
            ui.add_space(12.0);
            ui.label(egui::RichText::new("Export").strong().size(13.0));
            ui.separator();

            ui.horizontal(|ui| {
                if ui.button("IBSimu CSV").clicked() {
                    let (beam, geo) = self.make_beam_and_geo();
                    match Tracker::export_to_ibsimu(&beam, &geo) {
                        std::result::Result::Ok(_) => {
                            self.status = "Exported beam_tracing.csv".into()
                        }
                        Err(e) => self.status = format!("Export error: {}", e),
                    }
                }
                if ui.button("FEMM Lookup").clicked() {
                    let (beam, geo) = self.make_beam_and_geo();
                    match Tracker::export_femm_lookup(&beam, &geo) {
                        std::result::Result::Ok(_) => {
                            self.status = "Exported FEMM-Lookup.csv".into()
                        }
                        Err(e) => self.status = format!("Export error: {}", e),
                    }
                }
            });
        }
    }

    fn draw_results(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        // ── Tab bar ───────────────────────────────────────────
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.active_tab, Tab::Results, "Results");
            ui.selectable_value(&mut self.active_tab, Tab::Fluxes, "Fluxes");
            ui.selectable_value(&mut self.active_tab, Tab::Energies, "Energies");
        });
        ui.separator();

        if self.tracker.is_none() {
            ui.add_space(20.0);
            ui.label(egui::RichText::new("No results yet.").color(egui::Color32::GRAY));
            return;
        }

        egui::ScrollArea::vertical().show(ui, |ui| match self.active_tab {
            Tab::Results => self.draw_tab_results(ui),
            Tab::Fluxes => self.draw_tab_fluxes(ui),
            Tab::Energies => self.draw_tab_energies(ui),
        });
    }

    fn draw_tab_results(&self, ui: &mut egui::Ui) {
        let t = self.tracker.as_ref().unwrap();
        let bore_m = self.r_gap_mm * 1e-3;
        let clip_x = t.max_env_x > bore_m * 0.95;
        let clip_y = t.max_env_y > bore_m * 0.95;
        let avg = (t.x_f + t.y_f) / 2.0;
        let imbal = if avg > 1e-12 {
            (t.x_f - t.y_f).abs() / avg * 100.0
        } else {
            0.0
        };

        // ── Optimizer outputs ─────────────────────────────────
        ui.label(egui::RichText::new("Optimizer").strong());
        egui::Grid::new("opt_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("MMF₁ (outer)");
                ui.label(format!("{:.1} A·t", self.mmf1.unwrap_or(0.0)));
                ui.end_row();

                ui.label("MMF₂ (inner)");
                ui.label(format!("{:.1} A·t", self.mmf2.unwrap_or(0.0)));
                ui.end_row();

                ui.label("g₁ (outer)");
                ui.label(format!("{:.3} T/m", self.g1.unwrap_or(0.0)));
                ui.end_row();

                ui.label("g₂ (inner)");
                ui.label(format!("{:.3} T/m", self.g2.unwrap_or(0.0)));
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Beam spot ─────────────────────────────────────────
        ui.label(egui::RichText::new("Focal spot").strong());
        egui::Grid::new("spot_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Spot x");
                ui.label(format!("{:.3} mm", t.x_f * 1000.0));
                ui.end_row();

                ui.label("Spot y");
                ui.label(format!("{:.3} mm", t.y_f * 1000.0));
                ui.end_row();

                ui.label("x/y imbalance");
                ui.colored_label(
                    if imbal > 10.0 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::WHITE
                    },
                    format!("{:.1} %", imbal),
                );
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Envelope ──────────────────────────────────────────
        ui.label(egui::RichText::new("Envelope").strong());
        egui::Grid::new("env_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Max env x");
                ui.colored_label(
                    if clip_x {
                        egui::Color32::RED
                    } else {
                        egui::Color32::WHITE
                    },
                    format!(
                        "{:.2} mm{}",
                        t.max_env_x * 1000.0,
                        if clip_x { "  ⚠ CLIPPED" } else { "" }
                    ),
                );
                ui.end_row();

                ui.label("Max env y");
                ui.colored_label(
                    if clip_y {
                        egui::Color32::RED
                    } else {
                        egui::Color32::WHITE
                    },
                    format!(
                        "{:.2} mm{}",
                        t.max_env_y * 1000.0,
                        if clip_y { "  ⚠ CLIPPED" } else { "" }
                    ),
                );
                ui.end_row();

                ui.label("Total length");
                ui.label(format!("{:.1} mm", t.total_length * 1000.0));
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Pole tip fields ───────────────────────────────────
        ui.label(egui::RichText::new("Pole-tip fields").strong());
        let geo = self.make_geo();
        let b1 = geo.solve_b_pole(self.mmf1.unwrap_or(0.0));
        let b2 = geo.solve_b_pole(self.mmf2.unwrap_or(0.0));

        egui::Grid::new("btip_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("B_tip outer");
                ui.colored_label(
                    if b1 > self.b_sat * 0.9 {
                        egui::Color32::RED
                    } else if b1 > self.b_sat * 0.7 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::WHITE
                    },
                    format!("{:.3} T", b1),
                );
                ui.end_row();

                ui.label("B_tip inner");
                ui.colored_label(
                    if b2 > self.b_sat * 0.9 {
                        egui::Color32::RED
                    } else if b2 > self.b_sat * 0.7 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::WHITE
                    },
                    format!("{:.3} T", b2),
                );
                ui.end_row();

                ui.label("B_sat limit");
                ui.label(format!("{:.2} T", self.b_sat));
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Beam physics ──────────────────────────────────────
        ui.label(egui::RichText::new("Beam physics").strong());
        egui::Grid::new("phys_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Bρ");
                ui.label(format!("{:.4} T·m", beam_rigidity(self.energy_mev)));
                ui.end_row();

                ui.label("x crossovers");
                ui.label(if t.x_xover.is_empty() {
                    "none".into()
                } else {
                    t.x_xover
                        .iter()
                        .map(|z| format!("{:.1} mm", z * 1000.0))
                        .collect::<Vec<_>>()
                        .join(", ")
                });
                ui.end_row();

                ui.label("y crossovers");
                ui.label(if t.y_xover.is_empty() {
                    "none".into()
                } else {
                    t.y_xover
                        .iter()
                        .map(|z| format!("{:.1} mm", z * 1000.0))
                        .collect::<Vec<_>>()
                        .join(", ")
                });
                ui.end_row();
            });
    }

    fn draw_tab_fluxes(&self, ui: &mut egui::Ui) {
        ui.label(egui::RichText::new("Outer quads (Q1 / Q3)").strong());
        egui::Grid::new("flux1_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("φ total");
                ui.label(format!("{:.3} μWb", self.phi_total1.unwrap_or(0.0) * 1e6));
                ui.end_row();

                ui.label("φ gap");
                ui.label(format!("{:.3} μWb", self.phi_gap1.unwrap_or(0.0) * 1e6));
                ui.end_row();

                ui.label("φ leakage");
                ui.label(format!("{:.3} μWb", self.phi_leak1.unwrap_or(0.0) * 1e6));
                ui.end_row();

                let total = self.phi_total1.unwrap_or(1.0);
                let leak = self.phi_leak1.unwrap_or(0.0);
                ui.label("Leakage fraction");
                ui.label(format!("{:.1} %", (leak / total).abs() * 100.0));
                ui.end_row();
            });

        ui.add_space(10.0);

        ui.label(egui::RichText::new("Inner quad (Q2)").strong());
        egui::Grid::new("flux2_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("φ total");
                ui.label(format!("{:.3} μWb", self.phi_total2.unwrap_or(0.0) * 1e6));
                ui.end_row();

                ui.label("φ gap");
                ui.label(format!("{:.3} μWb", self.phi_gap2.unwrap_or(0.0) * 1e6));
                ui.end_row();

                ui.label("φ leakage");
                ui.label(format!("{:.3} μWb", self.phi_leak2.unwrap_or(0.0) * 1e6));
                ui.end_row();

                let total = self.phi_total2.unwrap_or(1.0);
                let leak = self.phi_leak2.unwrap_or(0.0);
                ui.label("Leakage fraction");
                ui.label(format!("{:.1} %", (leak / total).abs() * 100.0));
                ui.end_row();
            });
    }

    fn draw_tab_energies(&self, ui: &mut egui::Ui) {
        ui.label(egui::RichText::new("Outer quads (Q1 / Q3)").strong());
        egui::Grid::new("energy1_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("E total");
                ui.label(format!("{:.4} J", self.e_total1.unwrap_or(0.0)));
                ui.end_row();

                ui.label("E gap");
                ui.label(format!("{:.4} J", self.e_gap1.unwrap_or(0.0)));
                ui.end_row();

                ui.label("E leakage");
                ui.label(format!("{:.4} J", self.e_leak1.unwrap_or(0.0)));
                ui.end_row();
            });

        ui.add_space(10.0);

        ui.label(egui::RichText::new("Inner quad (Q2)").strong());
        egui::Grid::new("energy2_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("E total");
                ui.label(format!("{:.4} J", self.e_total2.unwrap_or(0.0)));
                ui.end_row();

                ui.label("E gap");
                ui.label(format!("{:.4} J", self.e_gap2.unwrap_or(0.0)));
                ui.end_row();

                ui.label("E leakage");
                ui.label(format!("{:.4} J", self.e_leak2.unwrap_or(0.0)));
                ui.end_row();
            });
    }

    fn draw_plot(&self, ui: &mut egui::Ui) {
        let Some(t) = &self.tracker else {
            ui.centered_and_justified(|ui| {
                ui.label(
                    egui::RichText::new("Set parameters and press Run")
                        .size(18.0)
                        .color(egui::Color32::GRAY),
                );
            });
            return;
        };

        let bore_mm = self.r_gap_mm;
        let total_mm = t.total_length * 1000.0;

        let x_pos: PlotPoints =
            t.z.iter()
                .zip(t.x.iter())
                .map(|(&z, &x)| [z * 1000.0, x * 1000.0])
                .collect();
        let x_neg: PlotPoints =
            t.z.iter()
                .zip(t.x.iter())
                .map(|(&z, &x)| [z * 1000.0, -x * 1000.0])
                .collect();
        let y_pos: PlotPoints =
            t.z.iter()
                .zip(t.y.iter())
                .map(|(&z, &y)| [z * 1000.0, y * 1000.0])
                .collect();
        let y_neg: PlotPoints =
            t.z.iter()
                .zip(t.y.iter())
                .map(|(&z, &y)| [z * 1000.0, -y * 1000.0])
                .collect();

        Plot::new("envelope")
            .legend(egui_plot::Legend::default())
            .x_axis_label("z (mm)")
            .y_axis_label("r (mm)")
            .show(ui, |plot_ui| {
                // ── Quad region boundaries ────────────────────
                let quad_col = egui::Color32::from_rgba_unmultiplied(80, 80, 220, 80);
                let quad_boundaries = [
                    (0.0, t.q1_end * 1000.0, "Q1"),
                    (t.q2_start * 1000.0, t.q2_end * 1000.0, "Q2"),
                    (t.q3_start * 1000.0, t.q3_end * 1000.0, "Q3"),
                ];

                for (z_start, z_end, label) in quad_boundaries {
                    // Left edge
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [z_start, -bore_mm * 1.2],
                            [z_start, bore_mm * 1.2],
                        ]))
                        .color(quad_col)
                        .width(1.5)
                        .name(label),
                    );

                    // Right edge
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [z_end, -bore_mm * 1.2],
                            [z_end, bore_mm * 1.2],
                        ]))
                        .color(quad_col)
                        .width(1.5),
                    );

                    // Top cap
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [z_start, bore_mm * 1.2],
                            [z_end, bore_mm * 1.2],
                        ]))
                        .color(quad_col)
                        .width(1.5),
                    );

                    // Bottom cap
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [z_start, -bore_mm * 1.2],
                            [z_end, -bore_mm * 1.2],
                        ]))
                        .color(quad_col)
                        .width(1.5),
                    );
                }

                // ── Bore limit ────────────────────────────────
                let bore_col = egui::Color32::from_rgb(200, 60, 60);
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![[0.0, bore_mm], [total_mm, bore_mm]]))
                        .color(bore_col)
                        .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                        .name(format!("Bore ±{:.1} mm", bore_mm))
                        .width(1.5),
                );

                plot_ui.line(
                    Line::new(PlotPoints::new(vec![[0.0, -bore_mm], [total_mm, -bore_mm]]))
                        .color(bore_col)
                        .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                        .width(1.5),
                );

                // ── Zero axis ─────────────────────────────────
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![[0.0, 0.0], [total_mm, 0.0]]))
                        .color(egui::Color32::from_rgba_unmultiplied(180, 180, 180, 60))
                        .width(1.0),
                );

                // ── Beam envelopes ────────────────────────────
                plot_ui.line(
                    Line::new(x_pos)
                        .color(egui::Color32::from_rgb(60, 130, 220))
                        .width(2.5)
                        .name("x-plane"),
                );
                plot_ui.line(
                    Line::new(x_neg)
                        .color(egui::Color32::from_rgb(60, 130, 220))
                        .width(2.5),
                );
                plot_ui.line(
                    Line::new(y_pos)
                        .color(egui::Color32::from_rgb(220, 80, 50))
                        .width(2.5)
                        .name("y-plane"),
                );
                plot_ui.line(
                    Line::new(y_neg)
                        .color(egui::Color32::from_rgb(220, 80, 50))
                        .width(2.5),
                );

                // ── Crossover markers ─────────────────────────
                for &zc in &t.x_xover {
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [zc * 1000.0, -bore_mm * 0.15],
                            [zc * 1000.0, bore_mm * 0.15],
                        ]))
                        .color(egui::Color32::from_rgb(60, 130, 220))
                        .width(1.5),
                    );
                }
                for &zc in &t.y_xover {
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![
                            [zc * 1000.0, -bore_mm * 0.15],
                            [zc * 1000.0, bore_mm * 0.15],
                        ]))
                        .color(egui::Color32::from_rgb(220, 80, 50))
                        .width(1.5),
                    );
                }
            });
    }
}

// ============================================================
// App logic
// ============================================================

impl QuadApp {
    fn make_geo(&self) -> MagnetGeometry {
        MagnetGeometry::new(
            self.r_gap_mm * 1e-3,
            self.l_mag_in * 0.0254,
            self.w_pole_mm * 1e-3,
            self.l_iron_mm * 1e-3,
            self.a_iron_mm2 * 1e-6,
            self.mu_i,
            self.b_sat,
            self.gap * 0.0254,
        )
    }

    fn make_beam_and_geo(&self) -> (Beam, MagnetGeometry) {
        let beam = Beam::new(
            self.drift_in * 0.0254,
            self.energy_mev,
            self.x0_mm * 1e-3,
            self.xp0_mrad * 1e-3,
        );
        (beam, self.make_geo())
    }

    fn run(&mut self) {
        self.status = "Running optimizer...".into();
        let (beam, geo) = self.make_beam_and_geo();

        match Tracker::optimize_mmf(&beam, &geo) {
            Some((mmf1, mmf2)) => {
                let g1 = geo.field_gradient(mmf1);
                let g2 = geo.field_gradient(mmf2);

                match Tracker::new(&beam, &geo, g1, g2, 400) {
                    std::result::Result::Ok(t) => {
                        // Flux results
                        let b1 = geo.solve_b_pole(mmf1);
                        let b2 = geo.solve_b_pole(mmf2);

                        let (phi1, phi_gap1, phi_leak1) = geo.solve_fluxes(b1, mmf1);
                        let (phi2, phi_gap2, phi_leak2) = geo.solve_fluxes(b2, mmf2);

                        self.phi_total1 = Some(phi1);
                        self.phi_gap1 = Some(phi_gap1);
                        self.phi_leak1 = Some(phi_leak1);
                        self.phi_total2 = Some(phi2);
                        self.phi_gap2 = Some(phi_gap2);
                        self.phi_leak2 = Some(phi_leak2);

                        // Energy results
                        let (e1, e_gap1, e_leak1) = geo.magnetic_energies(b1, mmf1);
                        let (e2, e_gap2, e_leak2) = geo.magnetic_energies(b2, mmf2);

                        self.e_total1 = Some(e1);
                        self.e_gap1 = Some(e_gap1);
                        self.e_leak1 = Some(e_leak1);
                        self.e_total2 = Some(e2);
                        self.e_gap2 = Some(e_gap2);
                        self.e_leak2 = Some(e_leak2);

                        self.mmf1 = Some(mmf1);
                        self.mmf2 = Some(mmf2);
                        self.g1 = Some(g1);
                        self.g2 = Some(g2);
                        self.tracker = Some(t);
                        self.status = "Done.".into();
                    }
                    Err(e) => self.status = format!("Tracker error: {}", e),
                }
            }
            None => self.status = "Optimizer did not converge.".into(),
        }
    }
}

// ============================================================
// Entry point
// ============================================================

pub fn launch_gui() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Quadrupole Triplet Design Tool")
            .with_inner_size([1400.0, 860.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Quadrupole Triplet Design Tool",
        options,
        Box::new(|_cc| Box::new(QuadApp::default())),
    )
}
