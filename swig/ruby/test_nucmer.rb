require 'minitest/autorun'
require 'mummer'

$bases = "ACGT"
def seq n
  s = ""
  n.times { s += $bases[rand(4)] }
  s
end

class TestNucmer < Minitest::Unit::TestCase
  def test_exact_pair
    s1 = seq(100)
    s2 = s1[80..-1] + seq(80)
    o = Mummer::Options.new.minmatch(10).mincluster(10)
    a = Mummer::align_sequences(s1, s2, o)

    assert_equal 1, a.size
    assert_equal 1, a[0].dirB
    assert_equal 81, a[0].sA
    assert_equal 100, a[0].eA
    assert_equal 1, a[0].sB
    assert_equal 20, a[0].eB
    assert_equal 0, a[0].Errors
    assert_equal 0, a[0].SimErrors
    assert_equal 0, a[0].delta.size
    assert_equal 1.0, a[0].identity
    assert_equal 1.0, a[0].similarity
  end
end
