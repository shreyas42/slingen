
For[ @it@; 0; @n@-(@nb@); @nb@ ]
{
	%% Equal( NL[ B_1 ], ldiv_lnn_ow( L, A_1 ) )
	@out0@[h(@m@,@m@,0), h(@nb@,@n@,@n@-(@it@+@nb@))] = ldiv_lnn_ow_opt(@m@, @nb@; @op0@[h(@m@,@m@,0), h(@m@,@m@,0)]#LowerTriangular#, @op1@[h(@m@,@m@,0), h(@nb@,@n@,@n@-(@it@+@nb@))]);
};

