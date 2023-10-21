use iced::Color;


#[derive(Debug, Clone)]
pub enum Cell {
    FreeCell,
    FluidCell {
        velocity: [f32; 2],
        pressure: f32,
        temperature: f32,
        rhs: f32,
        f: f32,
        g: f32
    },
    BoundaryConditionCell,
}

impl Cell {
    pub fn tick(&mut self, _: usize) {
        match self {
            Cell::FluidCell { pressure, .. } => {
                let to_add = 1.00 as f32;
                *pressure = *pressure + to_add;
            }
            _ => {}
        }
    }
}

pub fn color(cell: &Cell) -> Color {
    match cell {
        Cell::FreeCell => Color::WHITE,
        Cell::FluidCell { pressure, .. } => {
            let hue: f32 = pressure % 360.0;
            let saturation = 1.0;
            let lightness = 0.5;
            
            let (r, g, b) = hsl_to_rgb(hue, saturation, lightness);

            Color::from_rgba(r, g, b, 1.0)
        },
        Cell::BoundaryConditionCell => Color::BLACK
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