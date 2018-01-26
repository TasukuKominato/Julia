@printf "Start at %s\n" now()

include("sdpa_file.jl")
include("sdpa02_tsk.jl")

input_file = ["arch0.dat-s","arch2.dat-s","arch4.dat-s","arch8.dat-s","control1.dat-s","control2.dat-s","control3.dat-s"]
output_file = ["output-arch0-procs.dat" "output-arch2-procs.dat" "output-arch4-procs.dat" "output-arch8-procs.dat" "output-control1-procs.dat" "output-control2-procs.dat" "output-control3-procs.dat";
				"output-arch0-workers.dat" "output-arch2-workers.dat" "output-arch4-workers.dat" "output-arch8-workers.dat" "output-control1-workers.dat" "output-control2-workers.dat" "output-control3-workers.dat"]

for i=1:7
	@printf "%s\n" now()
	@printf "read %s\n" input_file[i]
	# read problem
	m, nBlock, blockStruct, b, F0, F = sdpa_read(input_file[i])
	
	# open procs output file
	fp = open(output_file[1,i], "w")
	@printf "open %s\n" output_file[1,i]

	# compute at nprocs = 1
	@printf "nprocs = %d\n" nprocs()
	v = sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
	all1 = v[3]
	schur1 = v[4]
	write(fp, "$(v[1]) $(v[2]) $(v[3]) $(v[4]) $(all1/v[3]) $(schur1/v[4])\n")

	# compute at nprocs = 2:32
	for p=0:4
		q=2^p
		addprocs(q)
		@printf "nprocs = %d\n" nprocs()
		v = sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
		write(fp, "$(v[1]) $(v[2]) $(v[3]) $(v[4]) $(all1/v[3]) $(schur1/v[4])\n")
	end

	# close procs output file
	close(fp)
	@printf "close %s\n" output_file[1,i]

	# reset nprocs to 1
	@printf "reset num of procs\n"
	N = nprocs()
	for j=2:N rmprocs(procs()[nprocs()]) end
	@printf "check: procs = %d\n" nprocs()

	# open workers output file
	fp = open(output_file[2,i], "w")
	@printf "open %s\n" output_file[2,i]

	# compute at nprocs = 1
	@printf "nprocs = %d\n" nprocs()
	v = sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
	all1 = v[3]
	schur1 = v[4]
	write(fp, "$(v[1]) $(v[2]) $(v[3]) $(v[4]) $(all1/v[3]) $(schur1/v[4])\n")

	# compute at nprocs = 3
	addprocs(2)
	@printf "nprocs = %d\n" nprocs()
	v = sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
	write(fp, "$(v[1]) $(v[2]) $(v[3]) $(v[4]) $(all1/v[3]) $(schur1/v[4])\n")

	# compute at nprocs = 5:33
	for p=1:4
		q=2^p
		addprocs(q)
		@printf "nprocs = %d\n" nprocs()
		v = sdpa_newton(m, nBlock, blockStruct, -b, -F0, -F)
		write(fp, "$(v[1]) $(v[2]) $(v[3]) $(v[4]) $(all1/v[3]) $(schur1/v[4])\n")
	end

	# close workers output file
	close(fp)
	@printf "close %s\n" output_file[2,i]

	# reset nprocs to 1
	@printf "reset num of procs\n"
	N = nprocs()
	for j=2:N rmprocs(procs()[nprocs()]) end
	@printf "check: procs = %d\n" nprocs()
end
