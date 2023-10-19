
use iced::Color;


#[derive(Clone, Debug)]
pub enum Cell {
    FreeCell,
    FluidCell {
        vertical_velocity: f32,
        horizontal_velocity: f32,
        pressure: f32,
        temperature: f32,
        rhs: f32,
        f: f32,
        g: f32
    },
    BoundaryConditionCell,
}

pub fn color(cell: &Cell) -> Color {
    match cell {
        Cell::FreeCell => Color::WHITE,
        Cell::FluidCell { pressure, .. } => {
            // Map the pressure to a full color spectrum
            let hue: f32 = (pressure).max(0.0).min(360.0); // Map pressure to hue in degrees (0-360)
            let saturation = 1.0; // Maximum saturation
            let lightness = 0.5; // Medium lightness (you can adjust this)

            // Convert HSL to RGB
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