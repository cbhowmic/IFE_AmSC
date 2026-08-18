SELECT
    a.name AS archive,
    d.datasetid,
    d.name AS run_id,
    CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) AS minimum_startup_inventory_g
FROM datasets d
JOIN archives a
    ON d.archiveid = a.archiveid
JOIN attributes at
    ON at.archiveid = d.archiveid
   AND at.datasetid = d.datasetid
WHERE a.name LIKE :archive_name
  AND at.name = '/output:I_startup (g)'
ORDER BY a.name, d.datasetid;
