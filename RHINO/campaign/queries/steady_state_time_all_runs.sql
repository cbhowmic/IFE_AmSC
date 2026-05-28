-- 
--  Show the time to ss for every run

SELECT
    a.name AS archive,
    d.name AS run_name,
    at.value AS output:Steady state time (days) 
FROM datasets d
JOIN archives a
    ON d.archiveid = a.archiveid
JOIN attributes at
    ON at.datasetid = d.datasetid
WHERE at.name = 'output:Steady state time (days)'
ORDER BY a.name, d.name;