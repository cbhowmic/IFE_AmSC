--
-- List runs with Tritium steady-state inventory and steady-state time

SELECT
  d.datasetid,
  d.dsid,
  d.name AS dataset,
  a.name AS archive,
  CAST(tri.value AS REAL) AS tritium_inventory,
  CAST(REPLACE(REPLACE(ss.value, '[', ''), ']', '') AS REAL) AS steady_state_time_days
FROM datasets d
JOIN archives a
  ON d.archiveid = a.archiveid
JOIN attributes tri
  ON tri.archiveid = d.archiveid
 AND tri.datasetid = d.datasetid
 AND tri.name = '/data/inventory/Tritium/mass_steady_total'
JOIN attributes ss
  ON ss.archiveid = d.archiveid
 AND ss.datasetid = d.datasetid
 AND ss.name = '/output:Steady state time (days)'
WHERE CAST(tri.value AS REAL) > 2
ORDER BY CAST(tri.value AS REAL) DESC;