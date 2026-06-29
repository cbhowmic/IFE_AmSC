-- 
--  Startup inventory greater than 1000 

SELECT
  d.datasetid,
  d.dsid,
  d.name AS dataset,
  a.name AS archive,
  CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) AS I_startup_g
FROM datasets d
JOIN archives a
  ON d.archiveid = a.archiveid
JOIN attributes at
  ON at.archiveid = d.archiveid
 AND at.datasetid = d.datasetid
WHERE a.name LIKE 'IFE/rhino%.aca'
  AND at.name = '/output:I_startup (g)'
  AND CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) > 2000
ORDER BY I_startup_g DESC;