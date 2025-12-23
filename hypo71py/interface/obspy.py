from pathlib import Path
from copy import deepcopy
from obspy import read_inventory

from hypo71py.model.station_phase import Station
from hypo71py.model.station_phase import build_pick_dict_from_event
from hypo71py.core.single import SINGLE

def load_stations_from_stationxml(path):
    """
    Load stations from StationXML.

    Parameters
    ----------
    path : str or Path
        Can be:
        - a single StationXML file
        - a directory containing StationXML files
        - a glob pattern

    Returns
    -------
    stations : list[Station]
    """
    path = Path(path)
    inventories = []

    if path.is_file():
        inventories.append(read_inventory(str(path)))

    elif path.is_dir():
        for f in sorted(path.glob("*.xml")):
            inventories.append(read_inventory(str(f)))

    else:
        # glob pattern
        for f in sorted(Path().glob(str(path))):
            inventories.append(read_inventory(str(f)))

    stations = []

    for inv in inventories:
        for net in inv:
            for sta in net:
                code = f"{net.code}.{sta.code}"
                lon = sta.longitude
                lat = sta.latitude
                elev = sta.elevation / 1000.0  # m → km
                stations.append(Station(code, lon, lat, elev))

    return stations

    


def event_to_hypo71_inputs(event):
    """
    Convert an ObsPy Event to HYPO71 inputs.

    Returns
    -------
    origin_info : dict
    pick_dict : dict
    picked_station_codes : set[str]
    """
    origin_info, pick_dict, _ = build_pick_dict_from_event(event)
    picked_station_codes = set(pick_dict.keys())
    return origin_info, pick_dict, picked_station_codes


def match_stations_to_picks(stations, picked_station_codes):
    """
    Filter stations to those used by picks.

    Returns
    -------
    stations_subset : list[Station]
    missing : set[str]
    """
    stations_subset = [
        s for s in stations
        if s.code.split('.')[-1] in picked_station_codes
    ]

    found = {s.code.split('.')[-1] for s in stations_subset}
    missing = picked_station_codes - found

    return stations_subset, missing


from obspy.core.event import Origin, OriginQuality, CreationInfo
from hypo71py.core.single import SINGLE


def relocate_event_obspy(
    event,
    stations,
    velocity_model,
    *,
    make_preferred=True,
    method_id="hypo71py",
    include_arrivals=False,
    use_fortran_speedups=False,
    verbose=False,
):
    origin_info, pick_dict, picked_codes = event_to_hypo71_inputs(event)
    stations_subset, missing = match_stations_to_picks(stations, picked_codes)

    loc = SINGLE(
        stations_subset,
        pick_dict,
        velocity_model,
        ZTR=origin_info["depth_km"],
        origin=(
            origin_info["longitude"],
            origin_info["latitude"],
            origin_info["time"],
        ),
        use_fortran_speedups=use_fortran_speedups,
        verbose=verbose,
    )

    lon, lat, depth, time, se, station_phases, qsd, ni, rms = loc

    nphases = int(station_phases.num_phases)
    nstations = int(station_phases.num_stations)

    origin = Origin(
        latitude=lat,
        longitude=lon,
        depth=depth * 1000.0,
        time=time,
        creation_info=CreationInfo(
            agency_id="hypo71py",
            method_id=method_id,
        ),
        quality=OriginQuality(
            standard_error=float(rms),
            used_phase_count=nphases,
            associated_phase_count=nphases,
            used_station_count=nstations,
        ),
)

    event.origins.append(origin)

    if make_preferred:
        event.preferred_origin_id = origin.resource_id

    return event, origin


    


def relocate_catalog_obspy(
    catalog,
    stations,
    velocity_model,
    *,
    make_preferred=True,
    method_id="hypo71py",
    include_arrivals=False,
    use_fortran_speedups=False,
    verbose=False,
    in_place=False,
    stop_on_error=False,
):
    """
    Relocate every Event in an ObsPy Catalog by attaching a new Origin.

    Parameters
    ----------
    catalog : obspy.core.event.Catalog
    stations : list[hypo71py.model.station_phase.Station]
    velocity_model : CrustalVelocityModel (or compatible)
    in_place : bool
        If True, modify the passed-in catalog. If False, work on a deepcopy.
    stop_on_error : bool
        If True, raise immediately on the first failing event.

    Returns
    -------
    cat_out : obspy.core.event.Catalog
        Catalog with new Origins appended to each successfully relocated event.
    summary : list[dict]
        Per-event summary (IDs, original vs relocated, rms, nphases, etc.)
    failures : list[dict]
        Per-event failure report (event_id, error string).
    """
    cat_out = catalog if in_place else deepcopy(catalog)

    summary = []
    failures = []

    for i, ev in enumerate(cat_out, 1):
        event_id = getattr(ev.resource_id, "id", None)

        # capture original preferred origin if present
        o0 = ev.preferred_origin() or (ev.origins[0] if ev.origins else None)
        orig_lat = getattr(o0, "latitude", None)
        orig_lon = getattr(o0, "longitude", None)
        orig_depth_km = (getattr(o0, "depth", None) / 1000.0) if getattr(o0, "depth", None) is not None else None
        orig_time = getattr(o0, "time", None)

        try:
            ev_out, new_origin = relocate_event_obspy(
                ev,
                stations,
                velocity_model,
                make_preferred=make_preferred,
                method_id=method_id,
                include_arrivals=include_arrivals,
                use_fortran_speedups=use_fortran_speedups,
                verbose=verbose,
            )

            # relocated values
            lat = new_origin.latitude
            lon = new_origin.longitude
            depth_km = (new_origin.depth / 1000.0) if new_origin.depth is not None else None
            t = new_origin.time

            q = new_origin.quality
            rms = getattr(q, "standard_error", None) if q is not None else None
            nph = getattr(q, "used_phase_count", None) if q is not None else None
            nst = getattr(q, "used_station_count", None) if q is not None else None
            gap = getattr(q, "azimuthal_gap", None) if q is not None else None

            summary.append({
                "i": i,
                "event_id": event_id,
                "orig_lat": orig_lat,
                "orig_lon": orig_lon,
                "orig_depth_km": orig_depth_km,
                "orig_time": orig_time,
                "reloc_lat": lat,
                "reloc_lon": lon,
                "reloc_depth_km": depth_km,
                "reloc_time": t,
                "rms": rms,
                "used_phase_count": nph,
                "used_station_count": nst,
                "azimuthal_gap": gap,
                "num_origins_after": len(ev_out.origins),
            })

        except Exception as e:
            failures.append({
                "i": i,
                "event_id": event_id,
                "error": f"{type(e).__name__}: {e}",
            })
            if stop_on_error:
                raise

    return cat_out, summary, failures