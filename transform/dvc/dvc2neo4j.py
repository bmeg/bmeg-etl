import yaml
import glob
import types
import os

""" simple transform of dvc into graph. Note: needs access to bolt://..."""

from neo4j import GraphDatabase

driver = GraphDatabase.driver(os.environ.get('NEO4J_URL', "bolt://localhost:7687"), auth=(os.environ['NEO4J_USERNAME'], os.environ['NEO4J_PASSWORD']))

path = '*.dvc'

with driver.session() as session:
    for filename in glob.iglob(path, recursive=True):
        with open(filename, 'r') as stream:
            dvc = yaml.load(stream)
            if 'cmd' not in dvc:
                dvc['cmd'] = ''
            dvc = types.SimpleNamespace(**dvc)
            # print(dvc.cmd)
            # print(dvc.md5)
            session.run(
                statement="MERGE (x:Command{md5: {md5} }) SET x = {dict_param}",
                parameters={'dict_param': {'cmd': dvc.cmd, 'md5': dvc.md5, 'filename': filename}, 'md5': dvc.md5}
            )
            if hasattr(dvc, 'deps'):
                for dep in dvc.deps:
                    print(filename, dep)
                    if 'md5' not in dep:
                        dep['md5'] = 'stub-md5.{}'.format(dep['path'])
                    session.run(
                        statement="MERGE (f:File{md5: {md5} }) SET f = {dict_param}",
                        parameters={'dict_param': dep, 'md5': dep['md5']}
                    )
                    session.run(
                        statement="MATCH (c:Command{md5: {cmd_md5} }) MATCH (f:File{md5: {file_md5} }) CREATE (c)-[:READS]->(f)",
                        parameters={'cmd_md5': dvc.md5, 'file_md5': dep['md5']}
                    )

            if hasattr(dvc, 'outs'):
                for out in dvc.outs:
                    print(filename, out)
                    if 'md5' not in out:
                        out['md5'] = 'stub-md5.{}'.format(out['path'])
                    session.run(
                        statement="MERGE (f:File{md5: {md5} }) SET f = {dict_param}",
                        parameters={'dict_param': out, 'md5': out['md5']}
                    )
                    session.run(
                        statement="MATCH (c:Command{md5: {cmd_md5} }) MATCH (f:File{md5: {file_md5} }) CREATE (c)-[:WRITES]->(f)",
                        parameters={'cmd_md5': dvc.md5, 'file_md5': out['md5']}
                    )
