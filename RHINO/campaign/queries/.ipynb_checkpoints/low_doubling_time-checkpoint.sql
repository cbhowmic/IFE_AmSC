--
-- Low doubling time

SELECT *
FROM (
  SELECT
    d.datasetid,
    d.dsid,
    d.name AS dataset,
    a.name AS archive,
    CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) AS doubling_time_days
  FROM datasets d
  JOIN archives a
    ON d.archiveid = a.archiveid
  JOIN attributes at
    ON at.archiveid = d.archiveid
   AND at.datasetid = d.datasetid
  WHERE a.name LIKE 'IFE/rhino%.aca'
    AND at.name = '/output:plant_doubling_time (days)'
)
WHERE doubling_time_days < 4
ORDER BY doubling_time_days ASC;