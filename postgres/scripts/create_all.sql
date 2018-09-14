CREATE TABLE vertex (
 gid varchar not null,
 label varchar not null,
 data json
);

CREATE TABLE edge (
 gid varchar not null,
 label varchar not null,
 from varchar not null,
 to varchar not null,
 data json
);
