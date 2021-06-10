class NovatelHeader:
    header_length = 0
    msg_name = ""
    msg_id = 0
    msg_type = 0
    port_add = 0
    port = ""
    msg_length = 0
    sequence = 0
    idle_time = 0
    time_status = 0
    gps_week = 0
    gps_sec = 0
    receiver_status = 0
    SW_version = 0


class NovatelBestpos:
    header = NovatelHeader()
    sol_status = 0
    pos_type = 0
    lat = 0
    lon = 0
    height = 0
    undulation = 0
    datum_id = 0
    lat_std = 0
    lon_std = 0
    height_std = 0
    base_id = 0
    diff_age = 0
    sol_age = 0
    sat_tracked_num = 0
    sat_used_num = 0
    sat_used_L1 = 0
    sat_used_L1L2 = 0
    reserved = 0
    ext_sol_status = 0
    gal_beidou_mask = 0
    gps_glo_mask = 0
