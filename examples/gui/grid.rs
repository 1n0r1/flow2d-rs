use flow2d_rs::cell::Cell;
use flow2d_rs::cell::CellType;
use flow2d_rs::simulation::Simulation;
use flow2d_rs::presets;

use iced::widget::canvas::{Cache, Canvas, Frame, Geometry, Path, Program, Stroke};
use iced::{mouse, Color, Element, Length, Point, Renderer, Size, Theme};

use plotters::prelude::*;

#[derive(Default)]
pub struct Grid {
    simulation: Simulation,
    next_cache: Cache,
    vector_cache: Cache,
    color_type: ColorType,
    zoom: f32,
    show_velocity: bool
}


#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Preset {
    #[default]
    CylinderCrossFlow,
    BackwardFacingStep,
    LidDrivenCavity
}

pub static ALLPRESET: &[Preset] = &[
    Preset::CylinderCrossFlow,
    Preset::BackwardFacingStep,
    Preset::LidDrivenCavity
];


#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum ColorType {
    #[default]
    Pressure,
    Speed,
    Streamline
}

pub static ALLCOLORTYPE: &[ColorType] = &[
    ColorType::Pressure,
    ColorType::Speed,
    ColorType::Streamline
];


impl Grid {
    pub fn set_preset(&mut self, preset: Preset) {
        self.next_cache.clear();
        self.vector_cache.clear();
        match preset {
            Preset::CylinderCrossFlow => {
                self.simulation = Simulation::from_preset(presets::cylinder_cross_flow());
            },
            Preset::BackwardFacingStep => {
                self.simulation = Simulation::from_preset(presets::backward_facing_step());
            },
            Preset::LidDrivenCavity => {
                self.simulation = Simulation::from_preset(presets::lid_driven_cavity());
            },
        }
        
    }
    
    pub fn set_zoom(&mut self, zoom: f32) {
        self.next_cache.clear();
        self.vector_cache.clear();
        self.zoom = zoom;
    }

    pub fn set_color_type(&mut self, color_type: ColorType) {
        self.next_cache.clear();
        self.vector_cache.clear();
        self.color_type = color_type
    }
    
    pub fn set_show_velocity(&mut self, show_velocity: bool) {
        self.next_cache.clear();
        self.vector_cache.clear();
        self.show_velocity = show_velocity;
    }

    pub fn get_time(&self) -> f32 {
        self.simulation.get_time()
    }

    pub fn view(&self) -> Element<()> {
        Canvas::new(self)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    pub fn tick(&mut self) {
        self.simulation.iterate_one_timestep();
        self.next_cache.clear();
        self.vector_cache.clear();
    }

    pub fn export_image(&self) {
        let scale = 20.0;
        let space_size = self.simulation.get_space_size();
        let delta_space = self.simulation.get_delta_space();

        // let pressure_range = self.simulation.get_pressure_range();
        let speed_range = self.simulation.get_speed_range();
        // let psi_range = self.simulation.get_psi_range();

        let pixel_scale = [scale * delta_space[0], scale * delta_space[1]];
        let drawing_area = BitMapBackend::new(
            "img/test.png",
            (
                ((space_size[0] as f32) * pixel_scale[0]) as u32,
                ((space_size[1] as f32) * pixel_scale[1]) as u32,
            ),
        )
        .into_drawing_area();
        drawing_area.fill(&WHITE).unwrap();

        for x in 0..self.simulation.get_space_size()[0] {
            for y in 0..self.simulation.get_space_size()[1] {
                let pos_x = x as i32;
                let reversed_y = self.simulation.get_space_size()[1] - 1 - y;
                let pos_y = reversed_y as i32;
                let color = color_speed(self.simulation.get_cell(x, y), speed_range);
                drawing_area
                    .draw(&Rectangle::new(
                        [
                            (
                                ((pos_x as f32) * pixel_scale[0]) as i32,
                                ((pos_y as f32) * pixel_scale[1]) as i32,
                            ),
                            (
                                ((pos_x as f32 + 1.0) * pixel_scale[0]) as i32,
                                ((pos_y as f32 + 1.0) * pixel_scale[1]) as i32,
                            ),
                        ],
                        Into::<ShapeStyle>::into(plotters::style::RGBColor(
                            (color.r * 255.0) as u8,
                            (color.g * 255.0) as u8,
                            (color.b * 255.0) as u8,
                        ))
                        .filled(),
                    ))
                    .unwrap();
            }
        }
        drawing_area.present().unwrap();
    }
}

impl Program<()> for Grid {
    type State = ();

    fn draw(
        &self,
        _state: &(),
        renderer: &Renderer,
        _theme: &Theme,
        bounds: iced::Rectangle,
        _cursor: mouse::Cursor,
    ) -> Vec<Geometry> {
        let cells = self.next_cache.draw(renderer, bounds.size(), |frame| {
            if self.zoom == 0.0 {
                frame.scale(100.0);
            } else {
                frame.scale(self.zoom);
            }

            self.draw_cells(frame);
        });
        
        let vectors = self.vector_cache.draw(renderer, bounds.size(), |frame| {
            if self.zoom == 0.0 {
                frame.scale(100.0);
            } else {
                frame.scale(self.zoom);
            }
            self.draw_velocity_vector(frame);
        });

        if self.show_velocity {
            vec![cells, vectors]
        } else {
            vec![cells]
        }
    }
}

impl Grid {
    fn draw_cells(&self, frame: &mut Frame) {
        let background = Path::rectangle(Point::ORIGIN, frame.size());
        frame.fill(&background, Color::from_rgb8(0x40, 0x44, 0x4B));

        frame.with_save(|frame| {
            let delta_x = self.simulation.get_delta_space()[0];
            let delta_y = self.simulation.get_delta_space()[1];
            let pressure_range = self.simulation.get_pressure_range();
            let speed_range = self.simulation.get_speed_range();
            let psi_range = self.simulation.get_psi_range();


            for x in 0..self.simulation.get_space_size()[0] {
                for y in 0..self.simulation.get_space_size()[1] {
                    let pos_x = delta_x * (x as f32);
                    let reversed_y = self.simulation.get_space_size()[1] - 1 - y;
                    let pos_y = delta_y * (reversed_y as f32);

                    let color: Color = match self.color_type {
                        ColorType::Pressure => color_presure(self.simulation.get_cell(x, y), pressure_range),
                        ColorType::Speed => color_speed(self.simulation.get_cell(x, y), speed_range),
                        ColorType::Streamline => color_psi(self.simulation.get_cell(x, y), psi_range),
                    };

                    frame.fill_rectangle(
                        Point::new(pos_x, pos_y),
                        Size::new(delta_x, delta_y),
                        color
                    );
                }
            }
        });
    }

    fn draw_velocity_vector(&self, frame: &mut Frame) {
        let delta_x = self.simulation.get_delta_space()[0];
        let delta_y = self.simulation.get_delta_space()[1];
        for x in 0..self.simulation.get_space_size()[0] {
            for y in 0..self.simulation.get_space_size()[1] {
                // if x % 2 != 0 || y % 2 != 0 {
                //     continue;
                // }
                if let CellType::FluidCell = self.simulation.get_cell(x, y).cell_type {
                    let pos_x = delta_x * (x as f32);
                    let reversed_y = self.simulation.get_space_size()[1] - 1 - y;
                    let pos_y = delta_y * (reversed_y as f32);

                    let velocity_scale = 0.1;
                    let velocity = self.simulation.get_centered_velocity(x, y);

                    let vector_start = Point::new(pos_x + delta_x / 2.0, pos_y + delta_y / 2.0);
                    let vector_end = Point::new(
                        vector_start.x + velocity[0] * velocity_scale,
                        vector_start.y - velocity[1] * velocity_scale,
                    );

                    let vector = Path::line(vector_start, vector_end);
                    frame.stroke(
                        &vector,
                        Stroke {
                            width: 1.0,
                            ..Default::default()
                        },
                    );
                }
            }
        }
    }
}

pub fn color_presure(cell: &Cell, pressure_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 = 240.0
                - (cell.pressure - pressure_range[0]) * 240.0
                    / (pressure_range[1] - pressure_range[0]);
            let saturation = 1.0;
            let lightness = 0.5;

            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)
        }
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}

pub fn color_speed(cell: &Cell, speed_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            let speed = (cell.velocity[0].powi(2) + cell.velocity[1].powi(2)).sqrt();

            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 =
                240.0 - (speed - speed_range[0]) * 240.0 / (speed_range[1] - speed_range[0]);
            let saturation = 1.0;
            let lightness = 0.5;

            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)
        }
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}

pub fn color_psi(cell: &Cell, psi_range: [f32; 2]) -> Color {
    match cell.cell_type {
        CellType::FluidCell => {
            // 240 offset to map from blue to red instead of the whole range of hue
            let hue: f32 =
                (cell.psi - psi_range[0]) * 1800.0 / (psi_range[1] - psi_range[0]) % 360.0;
            let saturation = 1.0;
            let lightness = 0.5;

            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgb(r, g, b)
        }
        CellType::BoundaryConditionCell(_) => Color::from_rgb(0.5, 0.5, 0.5),
        CellType::VoidCell => Color::BLACK,
    }
}

fn hsl_to_rgb(hue: f32, saturation: f32, lightness: f32) -> (f32, f32, f32) {
    let c = (1.0 - (2.0 * lightness - 1.0).abs()) * saturation;
    let x = c * (1.0 - ((hue / 60.0) % 2.0 - 1.0).abs());
    let m = lightness - c / 2.0;

    let (r, g, b) = if hue < 60.0 {
        (c, x, 0.0)
    } else if hue < 120.0 {
        (x, c, 0.0)
    } else if hue < 180.0 {
        (0.0, c, x)
    } else if hue < 240.0 {
        (0.0, x, c)
    } else if hue < 300.0 {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };

    (r + m, g + m, b + m)
}
