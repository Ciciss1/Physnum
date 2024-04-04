for i in range(nsimul_ntols):
    output_file = outputs_tols[i]
    cmd = f"{repertoire}{executable} {input_filename} adapt=true tol={tols[i]:.15g} output={output_file}"

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')