! UTEP Electronic Structure Lab (2020)
        integer :: group
        group = 111
        open(2,file='RUNS')
        write(2,*) 0,group
        write(2,*) 3,4
        write(2,*) 0
        close(2)

        end

