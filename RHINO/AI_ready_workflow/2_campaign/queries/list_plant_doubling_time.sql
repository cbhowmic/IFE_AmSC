SELECT
    a.name AS archive,
    d.datasetid,
    d.name AS run_id,
    CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) AS plant_doubling_time_days
FROM datasets d
JOIN archives a
    ON d.archiveid = a.archiveid
JOIN attributes at
    ON at.archiveid = d.archiveid
   AND at.datasetid = d.datasetid
WHERE at.name = '/output:plant_doubling_time (days)'
ORDER BY a.name, d.datasetid;