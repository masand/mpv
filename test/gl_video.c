#include "test_helpers.h"
#include "video/out/gl_video.h"

static void test_get_ambient_gamma(void **state) {
    float x;
    x = gl_video_scale_ambient_lux(16.0, 64.0, 2.40, 1.961, 16.0);
    assert_true(x == 2.40f);

    x = gl_video_scale_ambient_lux(16.0, 64.0, 2.40, 1.961, 64.0);
    assert_true(x == 1.961f);

    x = gl_video_scale_ambient_lux(16.0, 64.0, 1.961, 2.40, 64.0);
    assert_true(x == 2.40f);

    x = gl_video_scale_ambient_lux(16.0, 64.0, 2.40, 1.961, 0.0);
    assert_true(x == 2.40f);
}

int main(void) {
    const UnitTest tests[] = {
        unit_test(test_get_ambient_gamma),
    };
    return run_tests(tests);
}

