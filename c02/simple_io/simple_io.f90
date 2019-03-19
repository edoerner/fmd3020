program simple_io
    !
    ! Program that reads and print a name to console
    !
    implicit none

    character *20 :: first_name

    print *, 'type in your first name.' 
    print *, 'up to 20 characters'

    read *, first_name
    print *, first_name

end program simple_io