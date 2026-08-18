SELECT
    a.name AS archive,
    d.datasetid,
    d.name AS run_id,
    CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL)
        AS tritium_burning_rate_g_per_day
FROM datasets d
JOIN archives a
    ON d.archiveid = a.archiveid
JOIN attributes at
    ON at.archiveid = d.archiveid
   AND at.datasetid = d.datasetid
WHERE a.name LIKE :archive_name
  AND at.name = '/input:Ndotminus:Tritium burned per day'
ORDER BY a.name, d.datasetid;
