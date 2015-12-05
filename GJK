module gjk2d

	use types, only: sp
	
	implicit none
	
	integer :: maxInteraction = 30

contains

	function cross (a, b)
	
		real (sp) :: cross
		real (sp), intent(in) :: a (:), b (:)

		cross = a(1)*b(2) - a(1)*b(1)

	end function cross

	function triple_product (a,b,c)
	
		real (sp), intent(in) :: a (:), b (:), c (:)
		real (sp) :: triple_product (2)
		
		triple_product = b*dot_product(c,a) - a*dot_product(c,b)
	
	end function triple_product

	function support (A,n,V)
		
		integer, intent(in) :: n
		integer :: i
		real (sp), intent(in) :: A (:,:), V (:)
		real (sp) :: r , t, support (2)
		
		r = 0.0
		t = dot_product(A(1:2,1),V)
		support = A(1:2,1)
		
		do i = 1, n
			r = dot_product(A(1:2,i),V)
			if (r > t) then
				t = r
				support = A(1:2,i)
			end if
		end do
		
	end function support

	function support_difference (P,n,Q,m,V)
			
		integer, intent(in) :: n, m
		real (sp), intent(in) :: P (:, :), Q (:,:), V(:)
		real (sp) :: support_difference (2)
		
		support_difference = support(P,n,V) - support(Q,m,-V)
			
	end function support_difference
	
	subroutine simplex2 (S,V,ss,col) !line
	
		real (sp), intent(inout) :: S(:,:), V(:)
		integer, intent(inout) :: ss
		real (sp) :: a(2), b(2), ab(2), ao(2)
		logical, intent(out) :: col
		
		a = S(:,2)
		b = S(:,1)
		ab = b - a
		ao = -a
		V = triple_product(ab,ao,ab)
		col = .false.
		
	end subroutine simplex2 !line
	
	subroutine simplex3 (S,V,ss,col) !triangle
		
		real (sp), intent(inout) :: S(:,:), V(:)
		integer, intent(inout) :: ss
		real (sp) :: a(2), b(2), c(2), ab(2), ac(2), ao(2), abc(2), acb(2)
		logical, intent(out) :: col
		
		a = S(:,3)
		b = S(:,2)
		c = S(:,1)
		ao= -a
		ab = b - a
		ac = c - a
		abc = triple_product(ac,ab,ab)
		acb = triple_product(ab,ac,ac)
		
		if (dot_product(abc,ao) > 0) then 
			ss = 2
			S(:,2) = a
			S(:,1) = b
			V = abc
			col = .false.
		else
			if (dot_product(acb,ao) > 0) then
				ss = 2
				S(:,2) = a
				S(:,1) = c
				V = acb
				col = .false.
			else
				col = .true.
			end if
		end if

	end subroutine simplex3 !triangle
	
	subroutine doSimplex (S,V,col,ss)
	
		real (sp), intent(inout) :: S(:,:), V(:)
		integer, intent(inout) :: ss
		logical, intent(inout) :: col
	
		if (ss == 2) then !line segment
			call simplex2(S,V,ss,col)
		else if (ss == 3) then !triangle
			call simplex3(S,V,ss,col)
		end if
		
	end subroutine doSimplex
	
	logical function collision_2D (P,n,Q,m)
	
		integer, intent(in) :: n, m
		real (sp), intent(in) :: P (:, :), Q (:,:)
		
		real (sp) :: V (2), S (2,3)
		integer :: ss, counter
		logical :: col
		
		S(:,:) = 0.0
		V = P(1:2,1) - Q(1:2,1)
		S(1:2,1) = support_difference(P,n,Q,m,V)
		V = -V
		
		ss = 1
		counter = 1
		col = .false.
		collision_2D = .false.
		
		do
			ss = ss + 1
			S(1:2,ss) = support_difference(P,n,Q,m,V)
			if (dot_product(S(1:2,ss),V) < 0) then
				collision_2D = .false.
				exit
			else
				call doSimplex(S,V,col,ss)
				if(col) then 
					collision_2D = .true.
					exit
				end if
			end if
			
			if (counter > maxInteraction) then
				collision_2D = .true.
				write(*,'(a)') '  WARNING: out of interactions'
				exit
			end if
			
			counter = counter + 1
		end do
		
	end function collision_2D

end module GJK2D
