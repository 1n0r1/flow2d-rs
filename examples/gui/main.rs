mod grid;

use grid::Grid;

use std::time::{Duration, Instant};

use iced::executor;
use iced::theme::{self, Theme};
use iced::time;
use iced::widget::{button, column, container, row, slider, text};
use iced::window;
use iced::{Alignment, Application, Command, Element, Length, Settings, Subscription};

pub fn main() -> iced::Result {
    EulerFluidSimulation::run(Settings {
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
struct EulerFluidSimulation {
    grid: Grid,
    is_playing: bool,
    speed: usize,
}

#[derive(Debug, Clone)]
enum Message {
    Tick(Instant),
    TogglePlayback,
    Next,
    SpeedChanged(f32),
    Export,
    None,
}

impl Application for EulerFluidSimulation {
    type Message = Message;
    type Theme = Theme;
    type Executor = executor::Default;
    type Flags = ();

    fn new(_flags: ()) -> (Self, Command<Message>) {
        (
            Self {
                speed: 2,
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
                if !self.is_playing {
                    self.speed = speed.round() as usize;
                }
            }
            Message::Export => {
                self.grid.export_image();
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
        let selected_speed = self.speed;
        let controls = view_controls(self.is_playing, selected_speed);

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

fn view_controls<'a>(is_playing: bool, speed: usize) -> Element<'a, Message> {
    let playback_controls = row![
        button(if is_playing { "Pause" } else { "Play" }).on_press(Message::TogglePlayback),
        button("Next")
            .on_press(Message::Next)
            .style(theme::Button::Secondary),
        button("Export Image")
            .on_press(Message::Export)
            .style(theme::Button::Secondary),
    ]
    .spacing(10);

    let speed_controls = row![
        slider(1.0..=150.0, speed as f32, Message::SpeedChanged),
        text(format!("x{speed}")).size(16),
    ]
    .width(Length::Fill)
    .align_items(Alignment::Center)
    .spacing(10);

    row![playback_controls, speed_controls]
        .padding(10)
        .spacing(20)
        .align_items(Alignment::Center)
        .into()
}