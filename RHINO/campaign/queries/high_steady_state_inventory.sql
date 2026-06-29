-- 
--  Steady state inventory > 500

SELECT *
FROM (
  SELECT
    d.datasetid,
    d.dsid,
    d.name AS dataset,
    a.name AS archive,
    CAST(REPLACE(REPLACE(at.value, '[', ''), ']', '') AS REAL) AS steady_state_inventory_g
  FROM datasets d
  JOIN archives a
    ON d.archiveid = a.archiveid
  JOIN attributes at
    ON at.archiveid = d.archiveid
   AND at.datasetid = d.datasetid
  WHERE a.name LIKE 'IFE/rhino%.aca'
    AND at.name = '/output:Iops (g)'
)
WHERE steady_state_inventory_g > 4000
ORDER BY steady_state_inventory_g ASC;