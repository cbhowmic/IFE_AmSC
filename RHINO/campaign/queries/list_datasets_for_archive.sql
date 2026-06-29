-- 
--  Check the datasets under a specific archive

SELECT d.datasetid, d.dsid, d.name
FROM datasets d
JOIN archives a ON d.archiveid = a.archiveid
WHERE a.name = 'IFE/rhino1.aca'
ORDER BY d.datasetid
LIMIT 50;