# Security policy

Report security issues privately to the repository maintainers rather than in a
public issue.

PymoSAICS executes the program selected by the local user. Only select a MOSAICS
binary obtained from a trusted source. The plugin does not download executables,
invoke a shell, or elevate privileges.

Parameter input files can direct MOSAICS to read or write arbitrary paths that
the user account can access. Review an input file before running it.
