require 'test/unit'
require 'mummer'

$bases = "ACGT"
def seq n
  s = ""
  n.times { s += $bases[rand(4)] }
  s
end

def comp x
  x.tr("ACGTacgt", "TGCAtgca")
end

def rev_comp s
  s.reverse.tr("ACGTacgt", "TGCAtgca")
end

def assert_good_alignment(a, s1, s2)
  a.each { |al|
    errors = 0
    i, j = al.sA, (al.dirB == 1 ? al.sB : (s2.size - al.sB - 1))
    ni = i - 1
    ni += al.delta[0].abs - 1 unless al.delta.empty?
    k = 0

    while(i < al.eA) do
      b2 = al.dirB == 1 ? s2[j] : comp(s2[j])
      if i == ni && k < al.delta.size
        while i == ni && k < al.delta.size do
          if al.delta[k] > 0
            i += 1
          else
            j += al.dirB
          end
          errors += 1
          k += 1
          oni, ni = ni, i
          ni += al.delta[k].abs - 1 if k < al.delta.size
        end
      else
        errors += 1 if s1[i] != b2
        i += 1
        j += al.dirB
      end
    end
    assert_equal(al.eA, i)
    assert_equal(al.dirB == 1 ? al.eB : (s2.size - al.eB - 1), j)
    assert_equal(al.delta.size, k)
    assert_equal(al.SimErrors, errors)
  }
end

class TestNucmer < Test::Unit::TestCase
  def test_exact_pair
    s1 = seq(100)
    s2 = s1[80..-1] + seq(80)
    o = Mummer::Options.new.minmatch(10).mincluster(10)
    a = Mummer::align_sequences(s1, s2, o)

    assert(1 <= a.size)
    al = a.find { |x| x.sA == 81 && x.eA == 100 }
    refute_nil(al)
    assert_equal 1, al.dirB
    assert_equal 81, al.sA
    assert_equal 100, al.eA
    assert_equal 1, al.sB
    assert_equal 20, al.eB
    assert_equal 0, al.Errors
    assert_equal 0, al.SimErrors
    assert_equal 0, al.delta.size
    assert_equal 1.0, al.identity
    assert_equal 1.0, al.similarity
  end

  def test_pair
    s1 = seq(1000)
    s2 = s1[900..-1] + seq(900)
    s1[950] = ""
    s2[75] = ""
    s2[25] = s2[25] == "A" ? "C" : "A"
    rs1 = rev_comp(s1)
    rs2 = rev_comp(s2)
    assert_equal 999, s1.size
    assert_equal 999, s2.size
    o = Mummer::Options.new.minmatch(10).mincluster(10)

    a = Mummer::align_sequences(s1, s2, o)
    assert(1 <= a.size)
    al = a.find { |x| x.sA == 901 && x.eA == 999 }
    refute_nil al
    assert_equal 901, al.sA
    assert_equal 999, al.eA
    assert_equal 1, al.sB
    assert_equal 99, al.eB
    assert_equal 1, al.dirB
    assert_equal 3, al.Errors
    assert_equal 3, al.SimErrors
    assert_equal 2, al.delta.size
    assert_good_alignment(a, s1, s2)

    a = Mummer::align_sequences(s1, rs2, o)
    assert(1 <= a.size)
    al = a.find { |x| x.sA == 901 && x.eA == 999 }
    refute_nil al
    assert_equal 901, al.sA
    assert_equal 999, al.eA
    assert_equal 1, al.sB
    assert_equal 99, al.eB
    assert_equal -1, al.dirB
    assert_equal 3, al.Errors
    assert_equal 3, al.SimErrors
    assert_equal 2, al.delta.size
    assert_good_alignment(a, s1, rs2)

    a = Mummer::align_sequences(rs1, s2, o)
    assert(1 <= a.size)
    al = a.find { |x| x.sA == 1 && x.eA == 99 }
    refute_nil al
    assert_equal 1, al.sA
    assert_equal 99, al.eA
    assert_equal 901, al.sB
    assert_equal 999, al.eB
    assert_equal -1, al.dirB
    assert_equal 3, al.Errors
    assert_equal 3, al.SimErrors
    assert_equal 2, al.delta.size
    assert_good_alignment(a, rs1, s2)

    a = Mummer::align_sequences(rs1, rs2, o)
    assert(1 <= a.size)
    al = a.find { |x| x.sA == 1 && x.eA == 99 }
    refute_nil al
    assert_equal 1, al.sA
    assert_equal 99, al.eA
    assert_equal 901, al.sB
    assert_equal 999, al.eB
    assert_equal 1, al.dirB
    assert_equal 3, al.Errors
    assert_equal 3, al.SimErrors
    assert_equal 2, al.delta.size
    assert_good_alignment(a, rs1, rs2)
  end

  def test_threads
    nt = Mummer::get_num_threads
    assert_operator nt, :>=, 1
    Mummer::set_num_threads(nt + 1)
    nnt = Mummer::get_num_threads
    assert(nnt == 1 || nnt == nt + 1)
  end
end
