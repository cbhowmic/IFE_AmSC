--
-- List steady-state Tritium inventory for all runs

SELECT
    d.datasetid,
    d.dsid,
    d.name AS dataset,
    a.name AS archive,
    CAST(tri.value AS REAL) AS tritium_inventory
FROM datasets d
JOIN archives a
    ON d.archiveid = a.archiveid
JOIN attributes tri
    ON tri.archiveid = d.archiveid
   AND tri.datasetid = d.datasetid
WHERE tri.name = '/data/inventory/Tritium/mass_steady_total'
ORDER BY CAST(tri.value AS REAL) DESC;
