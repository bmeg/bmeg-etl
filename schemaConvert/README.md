## ETL scripts are run in this order from bmeg-etl root dir:

1. python schemaConvert/normalizeJSONSchema.py, using default arguments
2. python schemaEditTools/editSchema.py using hardcoded arugments

Note: for step 2 you will need the iceberg at schemaConvert/iceberg. You can get that on the main branch of:
    https://github.com/bmeg/iceberg


## Script Flow
1. Read from schema/ directory and write to schemaConvert/revisedSchemas
2. Read from schemaConvert/revisedSchemas and schemaConvert/iceberg and write to schemaConvert/schema_lifted
3. Old nodes like file, case, project, aliquot are removed from schemaConvert/revisedSchemas
4. The rest of the files form iceverg and revised schemas are moved to Destination directory
5. Files from schemaConvert/schema_lifted replace their older files in destination directory
6. schemaConvert/schema_lifted directory is removfed



