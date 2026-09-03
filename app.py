"""Solar Energy Planner - run with: python -m streamlit run app.py"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import requests
import streamlit as st
from scipy.optimize import differential_evolution


st.set_page_config(page_title="Solar Yield Optimizer", page_icon="☀️", layout="wide")
st.markdown("""
<style>
  .stApp {background: radial-gradient(circle at 78% 0%, #fa9b54 0%, #ec5d4e 20%, #45205f 52%, #101828 100%);}
  .block-container {padding-top: 2rem; max-width: 1280px;}
  h1, h2, h3, p, label {font-weight: 700 !important;}
  [data-testid="stSidebar"] {background: rgba(13, 18, 36, .93); border-right: 1px solid rgba(255, 255, 255, .18);}
  [data-testid="stSidebar"] * {color: #f8fafc !important;}
  div[data-testid="stMetric"] {background: rgba(255, 249, 237, .96); border: 1px solid #ffd27a; border-radius: 16px; padding: 18px; box-shadow: 0 8px 24px rgba(15, 23, 42, .22);}
  div[data-testid="stMetricLabel"] {font-size: .92rem; font-weight: 800; color: #334155 !important;}
  div[data-testid="stMetricValue"] {font-size: 2rem; font-weight: 900; color: #111827 !important;}
  [data-testid="stSidebar"] .stButton > button {background: linear-gradient(90deg, #f97316, #ef4444); color: black !important; border: 0; border-radius: 10px; font-weight: 900; min-height: 48px;}
  .hero {background: linear-gradient(105deg, rgba(15,23,42,.90), rgba(15,23,42,.46)); border: 1px solid rgba(255,255,255,.32); border-radius: 22px; padding: 26px 30px; margin: 0 0 22px 0; box-shadow: 0 12px 30px rgba(0,0,0,.22);}
  .hero h1 {color: #fff7ed; font-size: 2.65rem; margin: 0;}
  .hero p {color: #fff7ed; font-size: 1.08rem; margin: 8px 0 0;}
  .section-title {color: #fff7ed; font-size: 1.45rem; font-weight: 900; margin: 26px 0 8px;}
  [data-testid="stCaptionContainer"] p {color: #fff7ed !important;}
  .stAlert {border-radius: 12px; font-weight: 700;}
</style>
""", unsafe_allow_html=True)

NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"
FORECAST_URL = "https://api.open-meteo.com/v1/forecast"
HEADERS = {"User-Agent": "SolarEnergyPlanner/1.0 (student project)"}
ALBEDO = 0.25
SYSTEM_EFFICIENCY = 0.20  # Assumed 20% panel conversion efficiency.


def location_query(country: str, state: str, city: str, postal_code: str) -> str:
    """Create a location string while allowing any of the fields to be blank."""
    return ", ".join(value.strip() for value in [city, state, postal_code, country] if value.strip())


@st.cache_data(ttl=3600, show_spinner=False)
def geocode(query: str) -> tuple[float, float, str]:
    response = requests.get(
        NOMINATIM_URL,
        params={"q": query, "format": "json", "limit": 1},
        headers=HEADERS,
        timeout=15,
    )
    response.raise_for_status()
    results = response.json()
    if not results:
        raise ValueError("We could not find that location. Check the city, state, ZIP, and country.")
    result = results[0]
    return float(result["lat"]), float(result["lon"]), result["display_name"]


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_weather(latitude: float, longitude: float) -> tuple[pd.DataFrame, str]:
    params = {
        "latitude": latitude,
        "longitude": longitude,
        "hourly": "direct_normal_irradiance,diffuse_radiation,shortwave_radiation,cloud_cover,temperature_2m",
        "forecast_days": 7,
        "timezone": "auto",
    }
    response = requests.get(FORECAST_URL, params=params, timeout=20)
    response.raise_for_status()
    payload = response.json()
    if "hourly" not in payload:
        raise ValueError("The weather provider did not return an hourly forecast.")
    weather = pd.DataFrame(payload["hourly"])
    weather["time"] = pd.to_datetime(weather["time"])
    return weather, payload.get("timezone_abbreviation", "local time")


def demo_weather() -> pd.DataFrame:
    """A labelled fallback so the interface remains useful during an API outage."""
    times = pd.date_range(pd.Timestamp.now().floor("h"), periods=168, freq="h")
    daylight = np.clip(np.sin(np.pi * (times.hour.to_numpy() - 6) / 12), 0, None)
    return pd.DataFrame({
        "time": times,
        "direct_normal_irradiance": 800 * daylight,
        "diffuse_radiation": 110 * daylight,
        "shortwave_radiation": 920 * daylight,
        "cloud_cover": np.full(168, 15.0),
        "temperature_2m": 25 + 8 * daylight,
    })


def solar_position(times: pd.Series, latitude: float, longitude: float) -> tuple[np.ndarray, np.ndarray]:
    """Approximate local solar altitude and bearing; bearing is 0=south, -east, +west."""
    day = times.dt.dayofyear.to_numpy()
    local_hour = times.dt.hour.to_numpy() + times.dt.minute.to_numpy() / 60
    declination = np.radians(23.45 * np.sin(np.radians(360 * (284 + day) / 365)))
    solar_hour = local_hour + longitude / 15 - np.round(longitude / 15)
    hour_angle = np.radians(15 * (solar_hour - 12))
    phi = np.radians(latitude)
    altitude = np.arcsin(np.sin(phi) * np.sin(declination) + np.cos(phi) * np.cos(declination) * np.cos(hour_angle))
    bearing = np.arctan2(np.sin(hour_angle), np.cos(hour_angle) * np.sin(phi) - np.tan(declination) * np.cos(phi))
    return altitude, bearing


def plane_irradiance(tilt: float, azimuth: float, weather: pd.DataFrame, latitude: float, longitude: float) -> np.ndarray:
    """Calculate weather-adjusted irradiance on a tilted panel using a Liu-Jordan-style model."""
    altitude, sun_bearing = solar_position(weather["time"], latitude, longitude)
    tilt_r = np.radians(tilt)
    incidence = np.sin(altitude) * np.cos(tilt_r) + np.cos(altitude) * np.sin(tilt_r) * np.cos(sun_bearing - np.radians(azimuth))
    beam_ratio = np.maximum(incidence, 0) / np.maximum(np.sin(altitude), 0.065)
    dni = weather["direct_normal_irradiance"].to_numpy()
    diffuse = weather["diffuse_radiation"].to_numpy()
    ground = np.maximum(weather["shortwave_radiation"].to_numpy() - dni * np.maximum(np.sin(altitude), 0) - diffuse, 0)
    plane = dni * beam_ratio + diffuse * (1 + np.cos(tilt_r)) / 2 + ground * ALBEDO * (1 - np.cos(tilt_r)) / 2
    cloud_factor = 1 - weather["cloud_cover"].to_numpy() / 100 * 0.12
    soiling_factor = 1 - 0.12 * np.exp(-0.18 * tilt)
    return np.maximum(plane * cloud_factor * soiling_factor, 0)


def optimize(weather: pd.DataFrame, latitude: float, longitude: float) -> tuple[float, float]:
    def objective(values: np.ndarray) -> float:
        return -float(plane_irradiance(values[0], values[1], weather, latitude, longitude).sum())

    result = differential_evolution(objective, bounds=[(10, 60), (-45, 45)], seed=7, polish=True)
    return float(result.x[0]), float(result.x[1])


def direction_text(azimuth: float) -> str:
    if abs(azimuth) < 1:
        return "Face due south"
    side = "west" if azimuth > 0 else "east"
    return f"Face {abs(azimuth):.0f}° {side} of south"


def energy_table(irradiance: np.ndarray, weather: pd.DataFrame, area: float) -> pd.DataFrame:
    result = weather[["time"]].copy()
    result["panel_irradiance_wm2"] = irradiance
    result["energy_kwh"] = irradiance * area * SYSTEM_EFFICIENCY / 1000
    result["date"] = result["time"].dt.date
    return result


st.markdown("""
<div class="hero">
  <h1>☀️ Make your solar work harder.</h1>
  <p>Tell us where your panels are, and see the best angle, best collection times, and energy outlook.</p>
</div>
""", unsafe_allow_html=True)

with st.sidebar:
    st.header("1. Where are your panels?")
    country = st.text_input("Country", "United States")
    state = st.text_input("State / province", "Arizona")
    city = st.text_input("City", "Phoenix")
    postal_code = st.text_input("ZIP / postal code", "")
    st.divider()
    st.header("2. Your solar system")
    panel_area = st.slider("Panel area (m²)", 1.0, 100.0, 15.0, 0.5)
    interval = st.selectbox("How often can you adjust the panels?", ["Daily", "Monthly", "Quarterly", "Biannual", "Annual"], help="Uses live seven-day forecast data. Longer plans are forecast-based projections, not historical simulations.")
    run = st.button("Optimize my solar", type="primary", use_container_width=True)
    st.caption("Uses a 20% panel-efficiency estimate for simple planning.")

if "run" not in st.session_state:
    st.session_state.run = True
if run:
    st.session_state.run = True

if st.session_state.run:
    query = location_query(country, state, city, postal_code)
    if not query:
        st.error("Please enter at least a country, city, state/province, or ZIP/postal code.")
        st.stop()

    try:
        with st.spinner("Checking the location and loading the seven-day solar forecast..."):
            latitude, longitude, matched_location = geocode(query)
            weather, timezone = fetch_weather(latitude, longitude)
        source = f"Live Open-Meteo forecast ({timezone})"
    except (requests.RequestException, ValueError, KeyError) as error:
        latitude, longitude = 33.4484, -112.0740
        matched_location, timezone = "Phoenix, Arizona (demo location)", "local time"
        weather = demo_weather()
        source = "Demo forecast - live weather was unavailable"
        st.warning(f"{source}: {error}")

    optimization_weather = weather.iloc[:24].copy() if interval == "Daily" else weather
    best_tilt, best_azimuth = optimize(optimization_weather, latitude, longitude)
    optimized_irradiance = plane_irradiance(best_tilt, best_azimuth, weather, latitude, longitude)
    latitude_tilt = float(np.clip(abs(latitude), 10, 60))
    baseline_irradiance = plane_irradiance(latitude_tilt, 0, weather, latitude, longitude)
    optimized = energy_table(optimized_irradiance, weather, panel_area)
    baseline = energy_table(baseline_irradiance, weather, panel_area)

    daily_energy = optimized.groupby("date", as_index=False)["energy_kwh"].sum()
    predicted_daily = float(daily_energy.iloc[0]["energy_kwh"])
    baseline_daily = float(baseline.groupby("date")["energy_kwh"].sum().iloc[0])
    gain = (predicted_daily / baseline_daily - 1) * 100 if baseline_daily else 0
    days_in_plan = {"Daily": 1, "Monthly": 30, "Quarterly": 91, "Biannual": 182, "Annual": 365}[interval]
    projected_total = predicted_daily if interval == "Daily" else float(daily_energy["energy_kwh"].mean() * days_in_plan)

    today = optimized.iloc[:24].copy()
    peak_row = today.loc[today["panel_irradiance_wm2"].idxmax()]
    prime = today[today["panel_irradiance_wm2"] >= peak_row["panel_irradiance_wm2"] * 0.75]
    if peak_row["panel_irradiance_wm2"] > 0 and not prime.empty:
        prime_start, prime_end = prime.iloc[0]["time"], prime.iloc[-1]["time"] + pd.Timedelta(hours=1)
        prime_window = f"{prime_start.strftime('%I:%M %p').lstrip('0')} - {prime_end.strftime('%I:%M %p').lstrip('0')}"
        peak_time = peak_row["time"].strftime("%I:%M %p").lstrip("0")
    else:
        prime_start = prime_end = None
        prime_window, peak_time = "No strong solar window forecast", "No daylight forecast"

    st.success(f"Location found: {matched_location} | {latitude:.4f}°, {longitude:.4f}° | {source}")
    if interval != "Daily":
        st.info(f"{interval} total is a {days_in_plan}-day projection based on the average of this seven-day forecast.")

    st.markdown('<div class="section-title">Your best solar plan</div>', unsafe_allow_html=True)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Best panel tilt", f"{best_tilt:.1f}°", "Angle up from horizontal")
    c2.metric("Best direction to face", direction_text(best_azimuth), "South is usually strongest")
    c3.metric("Energy you could make today", f"{predicted_daily:.1f} kWh")
    c4.metric("Extra energy vs. latitude tilt", f"{gain:+.1f}%")
    st.metric(f"Estimated {interval.lower()} energy", f"{projected_total:.1f} kWh", f"For your {panel_area:.1f} m² panel area")

    st.markdown('<div class="section-title">Best time to collect solar energy</div>', unsafe_allow_html=True)
    timer_left, timer_right = st.columns(2)
    timer_left.metric("Best collection window", prime_window)
    timer_right.metric("Peak solar time", peak_time, f"{peak_row['panel_irradiance_wm2']:.0f} W/m² at the panel")
    st.caption("The best collection window includes hours forecast to reach at least 75% of the day’s peak solar input.")

    chart_left, chart_right = st.columns(2)
    hourly = pd.DataFrame({"Time": weather.iloc[:24]["time"], "Optimized": optimized_irradiance[:24], "Latitude tilt": baseline_irradiance[:24]}).melt("Time", var_name="Setup", value_name="Irradiance (W/m²)")
    with chart_left:
        st.markdown('<div class="section-title">Today’s solar capture</div>', unsafe_allow_html=True)
        line = px.line(hourly, x="Time", y="Irradiance (W/m²)", color="Setup", template="plotly_white")
        if prime_start is not None:
            line.add_vrect(x0=prime_start, x1=prime_end, fillcolor="#fbbf24", opacity=0.2, line_width=0, annotation_text="Best collection window")
        st.plotly_chart(line, use_container_width=True)
    with chart_right:
        st.markdown('<div class="section-title">7-day energy forecast</div>', unsafe_allow_html=True)
        daily_energy["date"] = daily_energy["date"].astype(str)
        bar = px.bar(daily_energy, x="date", y="energy_kwh", labels={"date": "Date", "energy_kwh": "Energy (kWh)"}, template="plotly_white", color_discrete_sequence=["#f59e0b"])
        st.plotly_chart(bar, use_container_width=True)

    st.markdown('<div class="section-title">Future energy prediction</div>', unsafe_allow_html=True)
    projection = pd.DataFrame({"Period": ["Daily", "Monthly", "Quarterly", "Biannual", "Annual"], "Estimated energy (kWh)": [predicted_daily, predicted_daily * 30, predicted_daily * 91, predicted_daily * 182, predicted_daily * 365]})
    st.plotly_chart(px.bar(projection, x="Period", y="Estimated energy (kWh)", template="plotly_white", color_discrete_sequence=["#0ea5e9"]), use_container_width=True)

    download = optimized.copy()
    download["recommended_tilt_deg"] = best_tilt
    download["recommended_azimuth_deg"] = best_azimuth
    download["location"] = matched_location
    st.download_button("Download hourly solar plan (CSV)", download.to_csv(index=False).encode("utf-8"), "solar_energy_plan.csv", "text/csv")
    st.caption("Planning estimate only. Verify roof structure, shading, electrical design, permits, and manufacturer requirements with a qualified solar professional.")
