mod grid;

use grid::{ColorType, Grid, Preset, ALLCOLORTYPE, ALLPRESET};

use std::time::{Duration, Instant};

use iced::executor;
use iced::theme::{self, Theme};
use iced::time;
use iced::widget::{button, checkbox, column, container, pick_list, row, slider, text};
use iced::window;
use iced::{Alignment, Application, Command, Element, Length, Settings, Subscription};

pub fn main() -> iced::Result {
    Flow2dGUI::run(Settings {
        antialiasing: true,
        window: window::Settings {
            size: (1600, 900),
            position: window::Position::Centered,
            ..window::Settings::default()
        },
        ..Settings::default()
    })
}

#[derive(Default)]
struct Flow2dGUI {
    grid: Grid,
    is_playing: bool,
    is_velocity_enabled: bool,
    speed: usize,
    preset: Preset,
    color_type: ColorType,
    zoom: f32,
}

#[derive(Debug, Clone)]
enum Message {
    Tick(Instant),
    TogglePlayback,
    Next,
    SpeedChanged(f32),
    Export,
    PresetPicked(Preset),
    ToggleVelocity(bool),
    ColorTypePicked(ColorType),
    ZoomChanged(f32),
    None,
}

impl std::fmt::Display for Preset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Preset::CylinderCrossFlow => "Cylinder Cross Flow",
                Preset::BackwardFacingStep => "Backward Facing Step",
                Preset::LidDrivenCavity => "Lid Driven Cavity",
            }
        )
    }
}

impl std::fmt::Display for ColorType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                ColorType::Pressure => "Pressure",
                ColorType::Speed => "Speed",
                ColorType::Streamline => "Streamline",
            }
        )
    }
}

impl Application for Flow2dGUI {
    type Message = Message;
    type Theme = Theme;
    type Executor = executor::Default;
    type Flags = ();

    fn new(_flags: ()) -> (Self, Command<Message>) {
        (
            Self {
                speed: 2,
                zoom: 100.0,
                ..Self::default()
            },
            Command::none(),
        )
    }

    fn title(&self) -> String {
        String::from("Euler Fluid Simulation")
    }

    fn update(&mut self, message: Message) -> Command<Message> {
        match message {
            Message::Tick(_) | Message::Next => {
                self.grid.tick();
            }
            Message::TogglePlayback => {
                self.is_playing = !self.is_playing;
            }
            Message::SpeedChanged(speed) => {
                self.speed = speed.round() as usize;
            }
            Message::Export => {
                self.grid.export_image();
            }
            Message::PresetPicked(preset) => {
                self.preset = preset;
                self.grid.set_preset(preset);
            }
            Message::ZoomChanged(zoom) => {
                self.grid.set_zoom(zoom);
                self.zoom = zoom;
            }
            Message::ColorTypePicked(color_type) => {
                self.color_type = color_type;
                self.grid.set_color_type(color_type);
            }
            Message::ToggleVelocity(is_velocity_enabled) => {
                self.is_velocity_enabled = is_velocity_enabled;
                self.grid.set_show_velocity(is_velocity_enabled);
            }
            Message::None => {}
        }
        Command::none()
    }

    fn subscription(&self) -> Subscription<Message> {
        if self.is_playing {
            time::every(Duration::from_millis(1000 / self.speed as u64)).map(Message::Tick)
        } else {
            Subscription::none()
        }
    }

    fn view(&self) -> Element<Message> {
        let controls = view_controls(
            self.is_playing,
            self.is_velocity_enabled,
            self.speed,
            self.preset,
            self.color_type,
            self.zoom,
        );

        let content = column![
            self.grid.view().map(move |_| { Message::None }),
            text(format!("time: {:.4}s", self.grid.get_time())).size(16),
            controls,
        ];

        container(content)
            .width(Length::Fill)
            .height(Length::Fill)
            .into()
    }

    fn theme(&self) -> Theme {
        Theme::Dark
    }
}

fn view_controls<'a>(
    is_playing: bool,
    is_velocity_enabled: bool,
    speed: usize,
    preset: Preset,
    color_type: ColorType,
    zoom: f32,
) -> Element<'a, Message> {
    let playback_controls = row![
        button(if is_playing { "Pause" } else { "Play" }).on_press(Message::TogglePlayback),
        button("Next")
            .on_press(Message::Next)
            .style(theme::Button::Secondary),
        button("Export Image")
            .on_press(Message::Export)
            .style(theme::Button::Secondary),
        checkbox("Velocity", is_velocity_enabled, Message::ToggleVelocity)
            .size(16)
            .spacing(5)
            .text_size(16),
        pick_list(ALLPRESET, Some(preset), Message::PresetPicked)
            .padding(8)
            .text_size(16),
        pick_list(ALLCOLORTYPE, Some(color_type), Message::ColorTypePicked)
            .padding(8)
            .text_size(16),
    ]
    .spacing(10);

    let speed_controls = row![
        slider(1.0..=300.0, speed as f32, Message::SpeedChanged),
        text(format!("Speed x{speed}")).size(16),
    ]
    .width(Length::Fill)
    .align_items(Alignment::Center)
    .spacing(10);

    let zoom_controls = row![
        slider(1.0..=1000.0, zoom as f32, Message::ZoomChanged),
        text(format!("Zoom x{zoom}")).size(16),
    ]
    .width(Length::Fill)
    .align_items(Alignment::Center)
    .spacing(10);

    row![playback_controls, column!(speed_controls, zoom_controls)]
        .padding(10)
        .spacing(20)
        .align_items(Alignment::Center)
        .into()
}
