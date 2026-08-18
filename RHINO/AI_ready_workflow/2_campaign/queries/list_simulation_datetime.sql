WITH run_datetimes AS (
    SELECT
        a.name AS archive,
        d.datasetid,
        d.name AS run_id,
        TRIM(at.value, '[]"') AS simulation_datetime
    FROM datasets d
    JOIN archives a
        ON d.archiveid = a.archiveid
    JOIN attributes at
        ON at.archiveid = d.archiveid
       AND at.datasetid = d.datasetid
    WHERE a.name LIKE :archive_name
      AND at.name = '/date'
)
SELECT
    archive,
    datasetid,
    run_id,
    SUBSTR(simulation_datetime, 1, 10) AS simulation_date,
    SUBSTR(simulation_datetime, 12, 8) AS simulation_time,
    simulation_datetime
FROM run_datetimes
ORDER BY archive, datasetid;
