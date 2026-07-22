# Security policy

Report security issues privately to the repository maintainers rather than in a
public issue.

PymoSAICS executes the bundled MOSAICS builds or a program selected by the local
user. Only select a custom binary obtained from a trusted source. The plugin
does not download executables, invoke a shell, or elevate privileges. RCSB PDB
fetching is limited to `https://files.rcsb.org`, validates four-character IDs,
and rejects responses larger than 20 MB.

Parameter input files can direct MOSAICS to read or write arbitrary paths that
the user account can access. Review an input file before running it.
