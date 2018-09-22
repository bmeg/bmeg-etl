# AgensGraph

## Install it:
  https://hub.docker.com/r/bitnine/agensgraph/

## Load it:
  See [loader.py](loader.py) and [config.yml](config.yml)

##  Issues:

  * volume mapping of data directory to host not working
  https://github.com/bitnine-oss/agensgraph-docker#usage-docker
  " not a database cluster"
    * resolution - don't start that way - e.g. do not map to host volume /home/agens/AgensGraph/data
    ```
    volumes:
      - ./:/src
      # - ./postgres/agensgraph/data:/home/agens/AgensGraph/data    
    ```

  * by default the postgres instance in the supplied docker file `does not` connect to external hostss
    * resolution - log onto docker image, edit
       vi /home/agens/AgensGraph/data/pg_hba.conf
       vi /home/agens/AgensGraph/data/postgresql.conf
       ```
       pg_hba.conf
       host all all 0.0.0.0/0 trused
       postgresql.conf
       listen_addresses='*'
       ```

  * cypher weirdness
  ```
  agens=# MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return n.values.ENSG00000000003 ;
   ensg00000000003
  -----------------

  (1 row)

  Time: 3.768 ms
  agens=# MATCH (n:Expression {gid: 'Expression:gtex:GTEX-1117F-0226-SM-5GZZ7'} ) return n.values['ENSG00000000003'] ;
         values
  --------------------
   46.040000915527344
  (1 row)
  ```
