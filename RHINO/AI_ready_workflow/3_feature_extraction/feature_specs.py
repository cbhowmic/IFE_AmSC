INPUT_SPECS = [
    {
        "key": "tritium_burning_rate",
        "role": "input",
        "source": "campaign_sql",
        "query": "list_tritium_burning_rate.sql",
        "column": "tritium_burning_rate_g_per_day",
    },
    {
        "key": "burn_fraction",
        "role": "input",
        "source": "campaign_sql",
        "query": "list_burn_fraction.sql",
        "column": "burn_fraction",
    },
]

OUTPUT_SPECS = [
    {
        "key": "plant_doubling_time_days",
        "role": "output",
        "source": "campaign_sql",
        "query": "list_plant_doubling_time.sql",
        "column": "plant_doubling_time_days",
    },
    {
        "key": "minimum_startup_inventory_g",
        "role": "output",
        "source": "campaign_sql",
        "query": "list_minimum_startup_inventory.sql",
        "column": "minimum_startup_inventory_g",
    },
    {
        "key": "tritium_in_isotope_separation",
        "role": "output",
        "source": "campaign_adios",
        "variable": "/data/inventory/Tritium/mass_steady",
        "subsystem": "Isotope_Seperation",
        "index": 10,
        "column": "tritium_in_isotope_separation",
    },
]

FEATURE_SPECS = INPUT_SPECS + OUTPUT_SPECS