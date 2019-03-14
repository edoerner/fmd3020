program simple_arithmetic
!
!   This program reads in three numbers and sums and averages them
!
!   Reference: Chivers I., Sleightholme J. (2018) Introduction to Programming. 
!   In: Introduction to Programming with Fortran. Springer, Cham
!
    implicit none
    real :: n1, n2, n3, average = 0.0, total = 0.0
    integer :: n = 3

    print *, 'Type in three numbers '
    print *, 'Separated by spaces or commas :'
    read *, n1, n2, n3

    total = n1 + n2 + n3
    average = total/n

    print *, 'Total sum of the numbers is ', total
    print *, 'Average of the numbers is ', average

end program simple_arithmetic