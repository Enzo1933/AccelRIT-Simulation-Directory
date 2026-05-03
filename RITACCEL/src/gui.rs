#![allow(non_snake_case)]
use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::{
    beam::Beam,
    einzel::EinzelGeometry,
    magnet::MagnetGeometry,
    tracker::{EinzelTracker, QuadTracker},
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
    r_gap_mm: f64,
    l_mag_in: f64,
    w_pole_mm: f64,
    l_iron_mm: f64,
    a_iron_mm2: f64,
    mu_i: f64,
    b_sat: f64,
    gap: f64,

    // ── Einzel Lens ───────────────────────────────────────────
    einzel_u_outer: f64,    // Outer electrode voltage (V)
    einzel_u_mid: f64,      // Middle electrode voltage (V)
    einzel_l_mid_mm: f64,   // Middle electrode length (mm)
    einzel_r_mm: f64,       // Cylinder radius (mm)
    einzel_start_z_mm: f64, // Tracking start (mm)
    einzel_end_z_mm: f64,   // Tracking end (mm)
    einzel_dz_mm: f64,      // Step size (mm)

    // ── Quad Results ──────────────────────────────────────────
    tracker: Option<QuadTracker>,
    mmf1: Option<f64>,
    mmf2: Option<f64>,
    g1: Option<f64>,
    g2: Option<f64>,

    phi_total1: Option<f64>,
    phi_gap1: Option<f64>,
    phi_leak1: Option<f64>,
    phi_total2: Option<f64>,
    phi_gap2: Option<f64>,
    phi_leak2: Option<f64>,

    e_total1: Option<f64>,
    e_gap1: Option<f64>,
    e_leak1: Option<f64>,
    e_total2: Option<f64>,
    e_gap2: Option<f64>,
    e_leak2: Option<f64>,

    // ── Einzel Results ────────────────────────────────────────
    einzel_tracker: Option<EinzelTracker>,

    status: String,

    // ── UI State ──────────────────────────────────────────────
    active_tab: Tab,
    active_plot: PlotView,
}

#[derive(PartialEq)]
enum Tab {
    Results,
    Fluxes,
    Energies,
    Einzel,
}

#[derive(PartialEq)]
enum PlotView {
    QuadEnvelope,
    EinzelProfile,
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

            einzel_u_outer: 0.0,
            einzel_u_mid: -5000.0,
            einzel_l_mid_mm: 50.0,
            einzel_r_mm: 20.0,
            einzel_start_z_mm: -150.0,
            einzel_end_z_mm: 150.0,
            einzel_dz_mm: 0.5,

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

            einzel_tracker: None,

            status: "Set parameters and press Run.".into(),
            active_tab: Tab::Results,
            active_plot: PlotView::QuadEnvelope,
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
            ui.heading("RITACCEL");
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

        // ── Central panel: plot with view switcher ────────────
        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(
                    &mut self.active_plot,
                    PlotView::QuadEnvelope,
                    "⬛ Quad Envelope",
                );
                ui.selectable_value(
                    &mut self.active_plot,
                    PlotView::EinzelProfile,
                    "⚡ Einzel Profile",
                );
            });
            ui.separator();

            match self.active_plot {
                PlotView::QuadEnvelope => self.draw_plot(ui),
                PlotView::EinzelProfile => self.draw_einzel_plot(ui),
            }
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
                ui.add(egui::Slider::new(&mut self.x0_mm, 0.5..=20000.0).step_by(0.1));
                ui.end_row();

                ui.label("x′₀ (mrad)");
                ui.add(egui::Slider::new(&mut self.xp0_mrad, 1.0..=20000.0).step_by(1.0));
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
                ui.add(egui::Slider::new(&mut self.r_gap_mm, 10.0..=15000.0).step_by(1.0));
                ui.end_row();

                ui.label("Magnet length (in)");
                ui.add(egui::Slider::new(&mut self.l_mag_in, 1.0..=12.0).step_by(0.1));
                ui.end_row();

                ui.label("Pole width (mm)");
                ui.add(egui::Slider::new(&mut self.w_pole_mm, 10.0..=2000.0).step_by(1.0));
                ui.end_row();

                ui.label("Iron path length (mm)");
                ui.add(egui::Slider::new(&mut self.l_iron_mm, 50.0..=10000.0).step_by(5.0));
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

                ui.label("Inter-magnet Gap (in)");
                ui.add(egui::Slider::new(&mut self.gap, 1.0..=100.0).step_by(0.1));
                ui.end_row();
            });

        ui.add_space(16.0);

        // ── Quad run button ───────────────────────────────────
        if ui
            .add_sized(
                [280.0, 40.0],
                egui::Button::new(egui::RichText::new("Run Quad Optimizer").size(15.0)),
            )
            .clicked()
        {
            self.run();
        }

        ui.add_space(6.0);

        // Status line
        let status_color = if self.status.starts_with("Error") || self.status.starts_with("No") {
            egui::Color32::RED
        } else if self.tracker.is_some() || self.einzel_tracker.is_some() {
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
                    match QuadTracker::export_to_ibsimu(&beam, &geo) {
                        std::result::Result::Ok(_) => {
                            self.status = "Exported beam_tracing.csv".into()
                        }
                        Err(e) => self.status = format!("Export error: {}", e),
                    }
                }
            });
        }

        // ── Einzel Lens section ───────────────────────────────
        ui.add_space(16.0);
        ui.label(egui::RichText::new("Einzel Lens").strong().size(13.0));
        ui.separator();

        egui::Grid::new("einzel_grid")
            .num_columns(2)
            .spacing([8.0, 6.0])
            .show(ui, |ui| {
                ui.label("U_outer (V)");
                ui.add(
                    egui::Slider::new(&mut self.einzel_u_outer, -20000.0..=20000.0)
                        .step_by(100.0),
                );
                ui.end_row();

                ui.label("U_mid (V)");
                ui.add(
                    egui::Slider::new(&mut self.einzel_u_mid, -20000.0..=20000.0)
                        .step_by(100.0),
                );
                ui.end_row();

                ui.label("L_mid (mm)");
                ui.add(egui::Slider::new(&mut self.einzel_l_mid_mm, 5.0..=200.0).step_by(1.0));
                ui.end_row();

                ui.label("R cylinder (mm)");
                ui.add(egui::Slider::new(&mut self.einzel_r_mm, 5.0..=100.0).step_by(0.5));
                ui.end_row();

                ui.label("Track start z (mm)");
                ui.add(
                    egui::Slider::new(&mut self.einzel_start_z_mm, -500.0..=0.0).step_by(5.0),
                );
                ui.end_row();

                ui.label("Track end z (mm)");
                ui.add(egui::Slider::new(&mut self.einzel_end_z_mm, 0.0..=500.0).step_by(5.0));
                ui.end_row();

                ui.label("Step size dz (mm)");
                ui.add(egui::Slider::new(&mut self.einzel_dz_mm, 0.1..=5.0).step_by(0.1));
                ui.end_row();
            });

        ui.add_space(12.0);

        if ui
            .add_sized(
                [280.0, 36.0],
                egui::Button::new(egui::RichText::new("Run Einzel Simulation").size(14.0)),
            )
            .clicked()
        {
            self.run_einzel();
        }
    }

    fn draw_results(&mut self, ui: &mut egui::Ui) {
        ui.add_space(8.0);

        // ── Tab bar ───────────────────────────────────────────
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.active_tab, Tab::Results, "Results");
            ui.selectable_value(&mut self.active_tab, Tab::Fluxes, "Fluxes");
            ui.selectable_value(&mut self.active_tab, Tab::Energies, "Energies");
            ui.selectable_value(&mut self.active_tab, Tab::Einzel, "Einzel");
        });
        ui.separator();

        // Quad-specific tabs require a completed quad run
        let needs_quad = matches!(self.active_tab, Tab::Results | Tab::Fluxes | Tab::Energies);
        if needs_quad && self.tracker.is_none() {
            ui.add_space(20.0);
            ui.label(
                egui::RichText::new("Run the optimizer first.").color(egui::Color32::GRAY),
            );
            return;
        }

        egui::ScrollArea::vertical().show(ui, |ui| match self.active_tab {
            Tab::Results => self.draw_tab_results(ui),
            Tab::Fluxes => self.draw_tab_fluxes(ui),
            Tab::Energies => self.draw_tab_energies(ui),
            Tab::Einzel => self.draw_tab_einzel(ui),
        });
    }

    // ── Quad results tabs (unchanged from original) ───────────

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
        ui.label(egui::RichText::new("Beam physics").strong());
        egui::Grid::new("phys_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
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

    // ── Einzel results tab ────────────────────────────────────

    fn draw_tab_einzel(&self, ui: &mut egui::Ui) {
        let Some(et) = &self.einzel_tracker else {
            ui.add_space(20.0);
            ui.label(
                egui::RichText::new("Run Einzel Simulation to see results.")
                    .color(egui::Color32::GRAY),
            );
            ui.add_space(6.0);
            ui.label(
                egui::RichText::new("Set parameters below and press\n\"Run Einzel Simulation\".")
                    .color(egui::Color32::DARK_GRAY)
                    .size(11.0),
            );
            return;
        };

        // ── Lens configuration ────────────────────────────────
        ui.label(egui::RichText::new("Lens Configuration").strong());
        egui::Grid::new("einzel_cfg_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("U_outer");
                ui.label(format!("{:.0} V", self.einzel_u_outer));
                ui.end_row();

                ui.label("U_mid");
                ui.label(format!("{:.0} V", self.einzel_u_mid));
                ui.end_row();

                let du = self.einzel_u_mid - self.einzel_u_outer;
                ui.label("ΔU (mid − outer)");
                ui.label(format!("{:.0} V", du));
                ui.end_row();

                ui.label("L_mid");
                ui.label(format!("{:.1} mm", self.einzel_l_mid_mm));
                ui.end_row();

                ui.label("R cylinder");
                ui.label(format!("{:.1} mm", self.einzel_r_mm));
                ui.end_row();

                let mode = if du < 0.0 {
                    "Focusing (ΔU < 0)"
                } else if du > 0.0 {
                    "Defocusing (ΔU > 0)"
                } else {
                    "No field (ΔU = 0)"
                };
                ui.label("Mode");
                ui.colored_label(
                    if du < 0.0 {
                        egui::Color32::from_rgb(60, 200, 100)
                    } else if du > 0.0 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::GRAY
                    },
                    mode,
                );
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Focal output ──────────────────────────────────────
        ui.label(egui::RichText::new("Focal Output").strong());
        egui::Grid::new("einzel_focal_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Final radius r_f");
                ui.label(format!("{:.4} mm", et.r_f * 1000.0));
                ui.end_row();

                ui.label("Final divergence r′_f");
                ui.label(format!("{:.4} mrad", et.r_prime_f * 1000.0));
                ui.end_row();

                // Thin-lens focal length estimate: f = -r_in / r'_out (parallel-entry beam)
                if et.r_prime_f.abs() > 1e-12 {
                    let f_mm = -(self.x0_mm) / et.r_prime_f;
                    ui.label("Est. focal length");
                    ui.colored_label(
                        if f_mm > 0.0 {
                            egui::Color32::from_rgb(60, 200, 100)
                        } else {
                            egui::Color32::YELLOW
                        },
                        format!("{:.1} mm", f_mm),
                    );
                    ui.end_row();
                }
            });

        ui.add_space(10.0);

        // ── Envelope health ───────────────────────────────────
        ui.label(egui::RichText::new("Beam Envelope").strong());
        let max_r_mm = et
            .r_phys
            .iter()
            .map(|r| r.abs())
            .fold(0.0_f64, f64::max)
            * 1000.0;
        let r_cyl = self.einzel_r_mm;
        let fill_pct = (max_r_mm / r_cyl) * 100.0;
        let clipped = max_r_mm > r_cyl;

        egui::Grid::new("einzel_env_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Max beam radius");
                ui.colored_label(
                    if clipped {
                        egui::Color32::RED
                    } else if fill_pct > 70.0 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::WHITE
                    },
                    format!("{:.3} mm", max_r_mm),
                );
                ui.end_row();

                ui.label("Cylinder fill");
                ui.colored_label(
                    if clipped {
                        egui::Color32::RED
                    } else if fill_pct > 70.0 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::from_rgb(60, 200, 100)
                    },
                    format!("{:.1} %{}", fill_pct, if clipped { "  ⚠ CLIPPED" } else { "" }),
                );
                ui.end_row();

                ui.label("Cylinder radius");
                ui.label(format!("{:.1} mm", r_cyl));
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── Voltage stats ─────────────────────────────────────
        ui.label(egui::RichText::new("On-axis Voltage").strong());
        let geo = self.make_einzel_geo();
        let peak_v = et
            .z
            .iter()
            .map(|&z| geo.voltage(z).abs())
            .fold(0.0_f64, f64::max);
        egui::Grid::new("einzel_v_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Applied U_mid");
                ui.label(format!("{:.0} V", self.einzel_u_mid));
                ui.end_row();

                ui.label("Peak |V(z)|");
                ui.label(format!("{:.1} V", peak_v));
                ui.end_row();

                // Show what fraction of beam energy the lens voltage is
                let beam_energy_ev = self.energy_mev * 1e6;
                let perturbation_pct = (peak_v / beam_energy_ev) * 100.0;
                ui.label("V / E_beam");
                ui.colored_label(
                    if perturbation_pct > 10.0 {
                        egui::Color32::YELLOW
                    } else {
                        egui::Color32::WHITE
                    },
                    format!("{:.3} %", perturbation_pct),
                );
                ui.end_row();
            });

        ui.add_space(10.0);

        // ── E-field stats ─────────────────────────────────────
        ui.label(egui::RichText::new("On-axis E-field").strong());

        // Find peak |E| and its z position
        let (peak_e, peak_e_z_mm) = et
            .e_field
            .iter()
            .zip(et.z.iter())
            .fold((0.0_f64, 0.0_f64), |(best_e, best_z), (&e, &z)| {
                if e.abs() > best_e { (e.abs(), z * 1000.0) } else { (best_e, best_z) }
            });

        // FWHM-style field width: fraction of z range where |E| > half peak
        let half_peak = peak_e * 0.5;
        let field_width_mm = {
            let active: Vec<f64> = et
                .z
                .iter()
                .zip(et.e_field.iter())
                .filter(|(_, e)| e.abs() >= half_peak)
                .map(|(&z, _)| z * 1000.0)
                .collect();
            if active.len() >= 2 {
                active.last().unwrap() - active.first().unwrap()
            } else {
                0.0
            }
        };

        egui::Grid::new("einzel_e_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                ui.label("Peak |E(z)|");
                ui.label(format!("{:.2} V/m", peak_e));
                ui.end_row();

                ui.label("Peak at z");
                ui.label(format!("{:.2} mm", peak_e_z_mm));
                ui.end_row();

                ui.label("FWHM field width");
                ui.label(format!("{:.2} mm", field_width_mm));
                ui.end_row();

                // E-field at the lens centre (z = 0)
                let e_at_centre = et
                    .z
                    .iter()
                    .zip(et.e_field.iter())
                    .min_by(|(za, _), (zb, _)| za.abs().partial_cmp(&zb.abs()).unwrap())
                    .map(|(_, &e)| e)
                    .unwrap_or(0.0);
                ui.label("E at z = 0");
                ui.label(format!("{:.2} V/m", e_at_centre));
                ui.end_row();
            });
    }

    // ── Quad envelope plot (original) ─────────────────────────

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

        let x_pos: PlotPoints = t
            .z
            .iter()
            .zip(t.x.iter())
            .map(|(&z, &x)| [z * 1000.0, x * 1000.0])
            .collect();
        let x_neg: PlotPoints = t
            .z
            .iter()
            .zip(t.x.iter())
            .map(|(&z, &x)| [z * 1000.0, -x * 1000.0])
            .collect();
        let y_pos: PlotPoints = t
            .z
            .iter()
            .zip(t.y.iter())
            .map(|(&z, &y)| [z * 1000.0, y * 1000.0])
            .collect();
        let y_neg: PlotPoints = t
            .z
            .iter()
            .zip(t.y.iter())
            .map(|(&z, &y)| [z * 1000.0, -y * 1000.0])
            .collect();

        Plot::new("envelope")
            .legend(egui_plot::Legend::default())
            .x_axis_label("z (mm)")
            .y_axis_label("r (mm)")
            .show(ui, |plot_ui| {
                let quad_col = egui::Color32::from_rgba_unmultiplied(80, 80, 220, 80);
                let quad_boundaries = [
                    (0.0, t.q1_end * 1000.0, "Q1"),
                    (t.q2_start * 1000.0, t.q2_end * 1000.0, "Q2"),
                    (t.q3_start * 1000.0, t.q3_end * 1000.0, "Q3"),
                ];

                for (z_start, z_end, label) in quad_boundaries {
                    let h = bore_mm * 1.2;
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![[z_start, -h], [z_start, h]]))
                            .color(quad_col)
                            .width(1.5)
                            .name(label),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![[z_end, -h], [z_end, h]]))
                            .color(quad_col)
                            .width(1.5),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![[z_start, h], [z_end, h]]))
                            .color(quad_col)
                            .width(1.5),
                    );
                    plot_ui.line(
                        Line::new(PlotPoints::new(vec![[z_start, -h], [z_end, -h]]))
                            .color(quad_col)
                            .width(1.5),
                    );
                }

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
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![[0.0, 0.0], [total_mm, 0.0]]))
                        .color(egui::Color32::from_rgba_unmultiplied(180, 180, 180, 60))
                        .width(1.0),
                );

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

    // ── Einzel profile plot ───────────────────────────────────

    fn draw_einzel_plot(&self, ui: &mut egui::Ui) {
        let Some(et) = &self.einzel_tracker else {
            ui.centered_and_justified(|ui| {
                ui.label(
                    egui::RichText::new("Run Einzel Simulation to view beam profile")
                        .size(18.0)
                        .color(egui::Color32::GRAY),
                );
            });
            return;
        };

        let geo = self.make_einzel_geo();

        // Compute max r and normalised voltage scale so V(z) fits within ±70% of the beam band
        let max_r_mm = et
            .r_phys
            .iter()
            .map(|r| r.abs())
            .fold(0.0_f64, f64::max)
            * 1000.0;
        let max_v = et
            .z
            .iter()
            .map(|&z| geo.voltage(z).abs())
            .fold(0.0_f64, f64::max);
        let v_scale_mm_per_v = if max_v > 1e-6 {
            max_r_mm * 0.7 / max_v
        } else {
            1.0
        };

        // Scale E-field so its peak also sits at ~50% of the beam band (separate from V scale)
        let max_e = et
            .e_field
            .iter()
            .map(|e| e.abs())
            .fold(0.0_f64, f64::max);
        let e_scale_mm_per_vm = if max_e > 1e-6 {
            max_r_mm * 0.5 / max_e
        } else {
            1.0
        };

        // Beam envelope (mirrored, in mm)
        let r_pos: PlotPoints = et
            .z
            .iter()
            .zip(et.r_phys.iter())
            .map(|(&z, &r)| [z * 1000.0, r.abs() * 1000.0])
            .collect();
        let r_neg: PlotPoints = et
            .z
            .iter()
            .zip(et.r_phys.iter())
            .map(|(&z, &r)| [z * 1000.0, -(r.abs() * 1000.0)])
            .collect();

        // On-axis voltage profile (scaled to mm for overlay)
        let v_line: PlotPoints = et
            .z
            .iter()
            .map(|&z| [z * 1000.0, geo.voltage(z) * v_scale_mm_per_v])
            .collect();

        // On-axis E-field profile (scaled to mm for overlay)
        let e_line: PlotPoints = et
            .z
            .iter()
            .zip(et.e_field.iter())
            .map(|(&z, &e)| [z * 1000.0, e * e_scale_mm_per_vm])
            .collect();

        let r_cyl_mm = self.einzel_r_mm;
        let z_start_mm = self.einzel_start_z_mm;
        let z_end_mm = self.einzel_end_z_mm;
        let half_l = self.einzel_l_mid_mm / 2.0;
        // Outer cylinders rendered with the same length as the middle electrode,
        // capped so they don't run past the tracking window.
        let outer_len = self.einzel_l_mid_mm;
        let outer_l_start = (z_start_mm).max(-half_l - outer_len);
        let outer_r_end   = (z_end_mm).min( half_l + outer_len);

        // Helper: draw a filled-wall cylinder cross-section (top + bottom wall at ±r, left + right end caps)
        let draw_cyl = |plot_ui: &mut egui_plot::PlotUi,
                        z0: f64, z1: f64,
                        col: egui::Color32,
                        label: Option<&str>| {
            let wall_h = r_cyl_mm;
            let wall_thick = r_cyl_mm * 0.18; // visual wall thickness
            // Top wall
            let top_outer = wall_h;
            let top_inner = wall_h - wall_thick;
            // Bottom wall (mirrored)
            let bot_inner = -top_inner;
            let bot_outer = -top_outer;

            // Top wall rectangle (4 sides)
            let top_pts = vec![
                [z0, top_inner], [z0, top_outer],
                [z1, top_outer], [z1, top_inner], [z0, top_inner],
            ];
            let bot_pts = vec![
                [z0, bot_inner], [z0, bot_outer],
                [z1, bot_outer], [z1, bot_inner], [z0, bot_inner],
            ];

            let mut top_line = Line::new(PlotPoints::new(top_pts)).color(col).width(1.5);
            if let Some(l) = label { top_line = top_line.name(l); }
            plot_ui.line(top_line);
            plot_ui.line(Line::new(PlotPoints::new(bot_pts)).color(col).width(1.5));
        };

        Plot::new("einzel_envelope")
            .legend(egui_plot::Legend::default())
            .x_axis_label("z (mm)")
            .y_axis_label("r (mm)")
            .show(ui, |plot_ui| {
                // ── Outer electrodes (blue-grey) ──────────────────
                let outer_col = egui::Color32::from_rgba_unmultiplied(120, 160, 210, 180);
                draw_cyl(plot_ui, outer_l_start, -half_l, outer_col, Some("Outer (U_outer)"));
                draw_cyl(plot_ui,  half_l, outer_r_end,  outer_col, None);

                // ── Middle electrode (amber) ──────────────────────
                let mid_col = egui::Color32::from_rgba_unmultiplied(220, 160, 40, 200);
                draw_cyl(plot_ui, -half_l, half_l, mid_col, Some("Middle (U_mid)"));

                // ── Cylinder aperture limit (red dashed) ──────────
                let cyl_col = egui::Color32::from_rgb(200, 60, 60);
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![
                        [z_start_mm, r_cyl_mm],
                        [z_end_mm, r_cyl_mm],
                    ]))
                    .color(cyl_col)
                    .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                    .name(format!("Aperture ±{:.1} mm", r_cyl_mm))
                    .width(1.5),
                );
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![
                        [z_start_mm, -r_cyl_mm],
                        [z_end_mm, -r_cyl_mm],
                    ]))
                    .color(cyl_col)
                    .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                    .width(1.5),
                );

                // ── Zero axis ─────────────────────────────────────
                plot_ui.line(
                    Line::new(PlotPoints::new(vec![[z_start_mm, 0.0], [z_end_mm, 0.0]]))
                        .color(egui::Color32::from_rgba_unmultiplied(180, 180, 180, 60))
                        .width(1.0),
                );

                // ── On-axis voltage profile (amber dashed) ────────
                plot_ui.line(
                    Line::new(v_line)
                        .color(egui::Color32::from_rgb(220, 170, 40))
                        .style(egui_plot::LineStyle::Dashed { length: 6.0 })
                        .width(1.5)
                        .name("V(z) [scaled]"),
                );

                // ── On-axis E-field profile (cyan dotted) ─────────
                plot_ui.line(
                    Line::new(e_line)
                        .color(egui::Color32::from_rgb(80, 210, 230))
                        .style(egui_plot::LineStyle::Dashed { length: 3.0 })
                        .width(1.5)
                        .name("E(z) [scaled]"),
                );

                // ── Beam envelope (green) ─────────────────────────
                let env_col = egui::Color32::from_rgb(60, 200, 100);
                plot_ui.line(
                    Line::new(r_pos)
                        .color(env_col)
                        .width(2.5)
                        .name("r(z)"),
                );
                plot_ui.line(Line::new(r_neg).color(env_col).width(2.5));
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

    fn make_einzel_geo(&self) -> EinzelGeometry {
        EinzelGeometry::new(
            self.einzel_u_outer,
            self.einzel_u_mid,
            self.einzel_l_mid_mm * 1e-3,
            self.einzel_r_mm * 1e-3,
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

        match QuadTracker::optimize_mmf(&beam, &geo) {
            Some((mmf1, mmf2)) => {
                let g1 = geo.field_gradient(mmf1);
                let g2 = geo.field_gradient(mmf2);

                match QuadTracker::new(&beam, &geo, g1, g2, 400) {
                    std::result::Result::Ok(t) => {
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
                        self.active_plot = PlotView::QuadEnvelope;
                        self.active_tab = Tab::Results;
                        self.status = "Done.".into();
                    }
                    Err(e) => self.status = format!("Tracker error: {}", e),
                }
            }
            None => self.status = "Optimizer did not converge.".into(),
        }
    }

    fn run_einzel(&mut self) {
        let beam = Beam::new(
            self.drift_in * 0.0254,
            self.energy_mev,
            self.x0_mm * 1e-3,
            self.xp0_mrad * 1e-3,
        );
        let geo = self.make_einzel_geo();

        let start_z = self.einzel_start_z_mm * 1e-3;
        let end_z = self.einzel_end_z_mm * 1e-3;
        let dz = self.einzel_dz_mm * 1e-3;

        if start_z >= end_z {
            self.status = "Error: start z must be less than end z.".into();
            return;
        }
        if dz <= 0.0 {
            self.status = "Error: step size must be positive.".into();
            return;
        }

        let tracker = EinzelTracker::new(&beam, &geo, start_z, end_z, dz);
        self.einzel_tracker = Some(tracker);
        self.active_plot = PlotView::EinzelProfile;
        self.active_tab = Tab::Einzel;
        self.status = "Einzel simulation complete.".into();
    }
}

// ============================================================
// Entry point
// ============================================================

pub fn launch_gui() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Optimization Engine")
            .with_inner_size([1400.0, 860.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Optimization Engine",
        options,
        Box::new(|_cc| Box::new(QuadApp::default())),
    )
}