use eframe::egui;
use egui_plot::{Line, Plot, PlotPoints};

use crate::physics::{Beam, Tracker, beam_rigidity, field_gradient};

// ============================================================
// App state
// ============================================================

pub struct QuadApp {
    // ── Beam parameters ──────────────────────────────────────
    l_mag_in:   f64,
    gap_in:     f64,
    drift_in:   f64,
    energy_mev: f64,
    x0_mm:      f64,
    xp0_mrad:   f64,

    // ── Magnet parameters ────────────────────────────────────
    bore_in:    f64,
    n1:         usize,
    n2:         usize,
    mu_r:       f64,
    sat:        f64,

    // ── Results ──────────────────────────────────────────────
    tracker:    Option<Tracker>,
    i1:         Option<f64>,
    i2:         Option<f64>,
    g1:         Option<f64>,
    g2:         Option<f64>,
    status:     String,
}

impl Default for QuadApp {
    fn default() -> Self {
        Self {
            l_mag_in:   4.0,
            gap_in:     4.0,
            drift_in:   18.0,
            energy_mev: 1.0,
            x0_mm:      5.0,
            xp0_mrad:   30.0,
            bore_in:    3.0,
            n1:         400,
            n2:         800,
            mu_r:       2000.0,
            sat:        1.5,
            tracker:    None,
            i1:         None,
            i2:         None,
            g1:         None,
            g2:         None,
            status:     "Press Run to simulate.".into(),
        }
    }
}

// ============================================================
// eframe App impl
// ============================================================

impl eframe::App for QuadApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // ── Top panel: title ─────────────────────────────────
        egui::TopBottomPanel::top("title").show(ctx, |ui| {
            ui.heading("Quadrupole Triplet Design Tool — 1 MeV Proton Beam");
        });

        // ── Left panel: controls ─────────────────────────────
        egui::SidePanel::left("controls").min_width(280.0).show(ctx, |ui| {
            ui.add_space(8.0);
            ui.label(egui::RichText::new("Beam parameters").strong());
            ui.separator();

            egui::Grid::new("beam_grid")
                .num_columns(2)
                .spacing([12.0, 6.0])
                .show(ui, |ui| {
                    ui.label("Magnet length (in)");
                    ui.add(egui::Slider::new(&mut self.l_mag_in, 2.0..=8.0).step_by(0.1));
                    ui.end_row();

                    ui.label("Gap (in)");
                    ui.add(egui::Slider::new(&mut self.gap_in, 1.0..=12.0).step_by(0.1));
                    ui.end_row();

                    ui.label("Drift to target (in)");
                    ui.add(egui::Slider::new(&mut self.drift_in, 6.0..=36.0).step_by(0.5));
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
            ui.label(egui::RichText::new("Magnet parameters").strong());
            ui.separator();

            egui::Grid::new("magnet_grid")
                .num_columns(2)
                .spacing([12.0, 6.0])
                .show(ui, |ui| {
                    ui.label("Bore radius (in)");
                    ui.add(egui::Slider::new(&mut self.bore_in, 0.5..=6.0).step_by(0.1));
                    ui.end_row();

                    ui.label("Outer turns N1");
                    ui.add(egui::Slider::new(&mut self.n1, 100..=1200).step_by(50.0));
                    ui.end_row();

                    ui.label("Inner turns N2");
                    ui.add(egui::Slider::new(&mut self.n2, 100..=1200).step_by(50.0));
                    ui.end_row();

                    ui.label("μᵣ (pole material)");
                    ui.add(egui::Slider::new(&mut self.mu_r, 100.0..=5000.0).step_by(100.0));
                    ui.end_row();

                    ui.label("Saturation B (T)");
                    ui.add(egui::Slider::new(&mut self.sat, 0.5..=2.5).step_by(0.05));
                    ui.end_row();
                });

            ui.add_space(16.0);

            // ── Run button ───────────────────────────────────
            if ui.add_sized(
                [240.0, 36.0],
                egui::Button::new(
                    egui::RichText::new("Run optimizer").size(15.0)
                ),
            ).clicked() {
                self.run();
            }

            ui.add_space(8.0);
            ui.label(egui::RichText::new(&self.status).color(
                if self.status.starts_with("Error") || self.status.starts_with("No") {
                    egui::Color32::RED
                } else if self.tracker.is_some() {
                    egui::Color32::GREEN
                } else {
                    egui::Color32::GRAY
                }
            ));

            // ── Info panel ───────────────────────────────────
            if let Some(t) = &self.tracker {
                ui.add_space(12.0);
                ui.label(egui::RichText::new("Results").strong());
                ui.separator();

                let bore_m  = self.bore_in * 0.0254;
                let clip_x  = t.max_env_x > bore_m * 0.95;
                let clip_y  = t.max_env_y > bore_m * 0.95;

                egui::Grid::new("results_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Gradient g1");
                        ui.label(format!("{:.3} T/m", self.g1.unwrap_or(0.0)));
                        ui.end_row();

                        ui.label("Gradient g2");
                        ui.label(format!("{:.3} T/m", self.g2.unwrap_or(0.0)));
                        ui.end_row();

                        ui.label("Spot x");
                        ui.label(format!("{:.3} mm", t.x_f * 1000.0));
                        ui.end_row();

                        ui.label("Spot y");
                        ui.label(format!("{:.3} mm", t.y_f * 1000.0));
                        ui.end_row();

                        let avg = (t.x_f + t.y_f) / 2.0;
                        let imbal = if avg > 1e-12 {
                            (t.x_f - t.y_f).abs() / avg * 100.0
                        } else { 0.0 };

                        ui.label("x/y imbalance");
                        ui.label(format!("{:.1} %", imbal));
                        ui.end_row();

                        ui.label("Max env x");
                        ui.colored_label(
                            if clip_x { egui::Color32::RED } else { egui::Color32::WHITE },
                            format!("{:.2} mm", t.max_env_x * 1000.0),
                        );
                        ui.end_row();

                        ui.label("Max env y");
                        ui.colored_label(
                            if clip_y { egui::Color32::RED } else { egui::Color32::WHITE },
                            format!("{:.2} mm", t.max_env_y * 1000.0),
                        );
                        ui.end_row();
                    });

                // ── Hardware outputs ─────────────────────────
                ui.add_space(12.0);
                ui.label(egui::RichText::new("Hardware outputs").strong());
                ui.separator();

                let r      = self.bore_in * 0.0254;
                let g1     = self.g1.unwrap_or(0.0);
                let g2     = self.g2.unwrap_or(0.0);
                let i1     = self.i1.unwrap_or(0.0);
                let i2     = self.i2.unwrap_or(0.0);
                let b_tip1 = g1 * r;
                let b_tip2 = g2 * r;
                let ni1    = b_tip1 * r / crate::MU0;
                let ni2    = b_tip2 * r / crate::MU0;

                egui::Grid::new("hw_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Outer pole-tip B");
                        ui.label(format!("{:.1} mT", b_tip1 * 1000.0));
                        ui.end_row();

                        ui.label("Inner pole-tip B");
                        ui.label(format!("{:.1} mT", b_tip2 * 1000.0));
                        ui.end_row();

                        ui.label("Outer A·turns");
                        ui.label(format!("{:.0} A·t", ni1));
                        ui.end_row();

                        ui.label("Inner A·turns");
                        ui.label(format!("{:.0} A·t", ni2));
                        ui.end_row();

                        ui.label("Outer current I1");
                        ui.label(format!("{:.3} A", i1));
                        ui.end_row();

                        ui.label("Inner current I2");
                        ui.label(format!("{:.3} A", i2));
                        ui.end_row();

                        ui.label("Outer Bρ");
                        ui.label(format!("{:.4} T·m",
                            beam_rigidity(self.energy_mev)));
                        ui.end_row();
                    });

                // ── Export buttons ───────────────────────────
                ui.add_space(12.0);
                if ui.button("Export IBSimu CSV").clicked() {
                    let beam = self.make_beam();
                    let _ = Tracker::export_to_ibsimu(
                        &beam, self.n1, self.n2, r, self.mu_r, self.sat, beam.L_mag_m, beam.gap_m
                    );
                    self.status = "Exported beam_tracing.csv".into();
                }
                if ui.button("Export FEMM lookup").clicked() {
                    let beam = self.make_beam();
                    let _ = Tracker::export_femm_lookup(
                        &beam, self.n1, self.n2, r, self.mu_r, self.sat, beam.L_mag_m, beam.gap_m
                    );
                    self.status = "Exported FEMM-Lookup.csv".into();
                }
            }
        });

        // ── Central panel: envelope plot ─────────────────────
        egui::CentralPanel::default().show(ctx, |ui| {
            if let Some(t) = &self.tracker {
                let bore_m = self.bore_in * 0.0254;

                // Convert to mm for display
                let z_mm: Vec<f64> = t.z.iter().map(|&z| z * 1000.0).collect();
                let x_mm: Vec<f64> = t.x.iter().map(|&x| x * 1000.0).collect();
                let y_mm: Vec<f64> = t.y.iter().map(|&y| y * 1000.0).collect();

                let x_pos: PlotPoints = z_mm.iter().zip(x_mm.iter())
                    .map(|(&z, &x)| [z, x]).collect();
                let x_neg: PlotPoints = z_mm.iter().zip(x_mm.iter())
                    .map(|(&z, &x)| [z, -x]).collect();
                let y_pos: PlotPoints = z_mm.iter().zip(y_mm.iter())
                    .map(|(&z, &y)| [z, y]).collect();
                let y_neg: PlotPoints = z_mm.iter().zip(y_mm.iter())
                    .map(|(&z, &y)| [z, -y]).collect();

                let bore_mm  = bore_m * 1000.0;
                let total_mm = t.total_length * 1000.0;

                Plot::new("envelope")
                    .legend(egui_plot::Legend::default())
                    .x_axis_label("z (mm)")
                    .y_axis_label("r (mm)")
                    .show(ui, |plot_ui| {

                        // Quad region shading via vertical lines
                        let q_color = egui::Color32::from_rgba_unmultiplied(80, 80, 220, 200);

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [0.0,              -bore_mm * 1.2],
                            [0.0,               bore_mm * 1.2],
                        ])).color(q_color).name("Q1"));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [t.q1_end * 1000.0,  -bore_mm * 1.2],
                            [t.q1_end * 1000.0,   bore_mm * 1.2],
                        ])).color(q_color));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [t.q2_start * 1000.0, -bore_mm * 1.2],
                            [t.q2_start * 1000.0,  bore_mm * 1.2],
                        ])).color(q_color).name("Q2"));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [t.q2_end * 1000.0,  -bore_mm * 1.2],
                            [t.q2_end * 1000.0,   bore_mm * 1.2],
                        ])).color(q_color));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [t.q3_start * 1000.0, -bore_mm * 1.2],
                            [t.q3_start * 1000.0,  bore_mm * 1.2],
                        ])).color(q_color).name("Q3"));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [t.q3_end * 1000.0,  -bore_mm * 1.2],
                            [t.q3_end * 1000.0,   bore_mm * 1.2],
                        ])).color(q_color));

                        // Bore limits
                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [0.0, bore_mm], [total_mm, bore_mm],
                        ])).color(egui::Color32::from_rgb(200, 60, 60))
                          .style(egui_plot::LineStyle::Dashed { length: 8.0 })
                          .name(format!("Bore ±{:.1} mm", bore_mm)));

                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [0.0, -bore_mm], [total_mm, -bore_mm],
                        ])).color(egui::Color32::from_rgb(200, 60, 60))
                          .style(egui_plot::LineStyle::Dashed { length: 8.0 }));

                        // Zero axis
                        plot_ui.line(Line::new(PlotPoints::new(vec![
                            [0.0, 0.0], [total_mm, 0.0],
                        ])).color(egui::Color32::from_rgba_unmultiplied(200,200,200,80)));

                        // Beam envelopes
                        plot_ui.line(Line::new(x_pos)
                            .color(egui::Color32::from_rgb(50, 120, 220))
                            .width(2.0)
                            .name("x-plane"));
                        plot_ui.line(Line::new(x_neg)
                            .color(egui::Color32::from_rgb(50, 120, 220))
                            .width(2.0));
                        plot_ui.line(Line::new(y_pos)
                            .color(egui::Color32::from_rgb(220, 80, 50))
                            .width(2.0)
                            .name("y-plane"));
                        plot_ui.line(Line::new(y_neg)
                            .color(egui::Color32::from_rgb(220, 80, 50))
                            .width(2.0));

                        // Crossover markers
                        for &zc in &t.x_xover {
                            plot_ui.line(Line::new(PlotPoints::new(vec![
                                [zc * 1000.0, -1.0],
                                [zc * 1000.0,  1.0],
                            ])).color(egui::Color32::from_rgb(50, 120, 220))
                              .width(1.0));
                        }
                        for &zc in &t.y_xover {
                            plot_ui.line(Line::new(PlotPoints::new(vec![
                                [zc * 1000.0, -1.0],
                                [zc * 1000.0,  1.0],
                            ])).color(egui::Color32::from_rgb(220, 80, 50))
                              .width(1.0));
                        }
                    });
            } else {
                ui.centered_and_justified(|ui| {
                    ui.label(egui::RichText::new("Set parameters and press Run")
                        .size(18.0)
                        .color(egui::Color32::GRAY));
                });
            }
        });
    }
}

// ============================================================
// App logic
// ============================================================

impl QuadApp {
    fn make_beam(&self) -> Beam {
        Beam::new(
            self.l_mag_in  * 0.0254,
            self.gap_in    * 0.0254,
            self.drift_in  * 0.0254,
            self.energy_mev,
            self.x0_mm     * 1e-3,
            self.xp0_mrad  * 1e-3,
        )
    }

    fn run(&mut self) {
        let beam = self.make_beam();
        let r    = self.bore_in * 0.0254;

        self.status = "Running optimizer...".into();

        match Tracker::optimize_nr(&beam, self.n1, self.n2, r, self.mu_r, self.sat, beam.L_mag_m, beam.gap_m) {
            Some((i1, i2)) => {
                let g1 = field_gradient(i1, self.n1, r, self.mu_r, self.sat, beam.L_mag_m, beam.gap_m);
                let g2 = field_gradient(i2, self.n2, r, self.mu_r, self.sat, beam.L_mag_m, beam.gap_m);

                match Tracker::new(&beam, g1, g2, 400) {
                    Ok(t) => {
                        self.i1      = Some(i1);
                        self.i2      = Some(i2);
                        self.g1      = Some(g1);
                        self.g2      = Some(g2);
                        self.tracker = Some(t);
                        self.status  = "Done.".into();
                    }
                    Err(e) => self.status = format!("Error: {}", e),
                }
            }
            None => self.status = "No solution found.".into(),
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
            .with_inner_size([1280.0, 800.0]),
        ..Default::default()
    };

    eframe::run_native(
        "Quadrupole Triplet Design Tool",
        options,
        Box::new(|_cc| Box::new(QuadApp::default())),
    )
}